use std::io::{BufReader, Read};

use leftover::{LeftoverReader, LeftoverView};
use stringbytes::StringBytes;

use crate::domain::{FastqParser, ParsingStep, ReadingState};
use crate::{Error, FastqRecord, ParseError};

/// Read four-line FASTQ records through a [`std::io::BufReader`].
///
/// Handles LF and CRLF and skips whitespace-only lines between records.
/// Wrapped sequences and qualities are not supported. A nonempty final quality
/// line may omit its newline; an empty quality line requires its own line ending.
/// All reading methods preserve the original record bytes and line endings.
///
/// [`read_view`](Self::read_view) borrows records, copying only records that
/// span buffer refills into a reusable buffer. [`read_owned`](Self::read_owned)
/// returns records that can outlive the reader.
/// Stop reading after a parse or I/O error. Input may resume with a new record
/// after EOF; the physical line count continues.
///
/// ```
/// use biorill_fastx::FastqReader;
///
/// let input = b"@r note\r\nAC\r\n+r\r\n12\r\n";
/// let mut reader = FastqReader::new(input.as_slice());
/// let record = reader.read_view()?.unwrap();
/// assert_eq!(record.id(), b"r");
/// assert_eq!(record.sequence(), b"AC");
/// assert_eq!(record.into_parts().0, input);
/// assert!(reader.read_view()?.is_none());
/// # Ok::<(), biorill_fastx::Error>(())
/// ```
pub struct FastqReader<R: Read> {
    reader: LeftoverReader<BufReader<R>>,
    parser: FastqParser,
    state: ReadingState,
}

impl<R: Read> FastqReader<R> {
    /// Create a FASTQ reader using [`std::io::BufReader`]'s default capacity.
    pub fn new(source: R) -> Self {
        Self {
            reader: LeftoverReader::new(BufReader::new(source)),
            parser: FastqParser::new(),
            state: ReadingState::new(),
        }
    }

    /// Create a FASTQ reader with a buffer of at least `capacity` bytes.
    /// Zero uses one byte so that reading can make progress. Records may be
    /// larger than the buffer; their bytes are collected across refills.
    pub fn with_capacity(capacity: usize, source: R) -> Self {
        Self {
            reader: LeftoverReader::new(BufReader::with_capacity(capacity.max(1), source)),
            parser: FastqParser::new(),
            state: ReadingState::new(),
        }
    }

    /// Number of physical lines processed, including skipped blank lines and
    /// a final line without a newline. Repeated reads at EOF process a retained
    /// whitespace-only tail again and count it each time.
    pub fn offset(&self) -> usize {
        self.state.offset
    }

    /// Read and validate the next record, borrowing its bytes until the next read.
    /// Records within the current buffer are borrowed directly. Records spanning
    /// refills are copied into reusable storage, which may allocate as it grows.
    /// Returns `None` at EOF when no record remains.
    ///
    /// # Errors
    /// Returns input errors, including interruptions, or a parse error for an
    /// invalid or incomplete record. Stop reading after an error.
    pub fn read_view(&mut self) -> Result<Option<FastqRecord<&[u8]>>, Error> {
        self.read_with(LeftoverView::take_view)
    }

    /// Read and validate the next record into a vector that can outlive the reader.
    /// Copies the record bytes into a new vector, including any bytes collected
    /// across refills. The reader keeps its reusable storage for later records.
    /// Returns `None` at EOF when no record remains.
    ///
    /// # Errors
    /// Returns the same input and parse errors as [`Self::read_view`].
    pub fn read_owned(&mut self) -> Result<Option<FastqRecord<Vec<u8>>>, Error> {
        self.read_with(LeftoverView::take_owned)
    }

    fn read_with<'a, T: StringBytes>(
        &'a mut self,
        take: impl FnOnce(LeftoverView<'a>, usize) -> T,
    ) -> Result<Option<FastqRecord<T>>, Error> {
        let state = &mut self.state;
        // The line start survives refills within this read only. Accumulating
        // preserves its position in leftover + chunk; a later read starts over.
        let mut start = 0;
        loop {
            let view = self.reader.fill_view()?;
            // Both constructors use a nonzero buffer capacity. An empty
            // filled chunk therefore means the source reported EOF.
            let eof = view.buffer().is_empty();
            match self
                .parser
                .parse(&view, state, &mut start, eof)
                .map_err(|source| Error::Parse {
                    line: state.offset,
                    source,
                })? {
                ParsingStep::NeedMore => {
                    // Keep every earlier record line and the unfinished tail.
                    // Accumulating preserves the parser's logical ranges.
                    view.accumulate();
                }
                ParsingStep::Skip { consumed } => view.consume(consumed),
                ParsingStep::Complete { parts, consumed } => {
                    // Take from the view that was parsed. Filling again here
                    // could read the source again after completing a record at EOF.
                    return FastqRecord::try_from((take(view, consumed), parts))
                        .map(Some)
                        .map_err(|source| Error::Parse {
                            line: state.offset,
                            source: ParseError::InvalidRecord { source },
                        });
                }
                // EOF leaves the retained bytes and parsing state in place.
                ParsingStep::Eof => return Ok(None),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::FastqReader;
    use crate::domain::RecordState;
    use crate::{Error, ParseError};

    type Result<T> = std::result::Result<T, Error>;

    fn create_reader(data: &str) -> FastqReader<Cursor<&[u8]>> {
        let reader = Cursor::new(data.as_bytes());
        FastqReader::new(reader)
    }

    #[test]
    fn test_read_valid_record() -> Result<()> {
        let fastq_data = "@seq1 description\nATGC\n+\n!!!!\n@seq2\nGCGT\n+\n$$$$\n";

        let mut reader = create_reader(fastq_data);

        let record = reader.read_owned()?.expect("Should have a record");

        // Check that the FASTQ fields match the expected values
        assert_eq!(record.id(), b"seq1");
        assert_eq!(record.description().unwrap(), b"description");
        assert_eq!(record.sequence(), b"ATGC");
        assert_eq!(record.separator(), b"");
        assert_eq!(record.quality(), b"!!!!");

        let record = reader.read_owned()?.expect("Should have a second record");
        assert_eq!(record.id(), b"seq2");
        assert!(record.description().is_none());
        assert!(reader.read_owned()?.is_none());

        Ok(())
    }

    #[test]
    fn test_header_description_boundaries() -> Result<()> {
        for separator in [" ", "\t"] {
            for description in ["", "x", "description", "MIRE{UMI:ACGT}", " text\tmore"] {
                for newline in ["\n", "\r\n"] {
                    let data = format!(
                        "@seq1{separator}{description}{newline}ACGT{newline}+{newline}IIII{newline}"
                    );
                    let mut reader = create_reader(&data);
                    let record = reader.read_owned()?.expect("Should have a record");
                    let expected = if description.is_empty() {
                        None
                    } else {
                        Some(description.as_bytes())
                    };
                    assert_eq!(record.id(), b"seq1");
                    assert_eq!(record.description(), expected, "input: {data:?}");
                    assert!(reader.read_owned()?.is_none());
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_invalid_header() -> Result<()> {
        let fastq_data = "seq1 description\nATGC\n+\n!!!!\n";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_owned();

        assert!(result.is_err()); // Expect error: Invalid header
        Ok(())
    }

    #[test]
    fn test_incomplete_record() -> Result<()> {
        let fastq_data = "@seq1 description\nATGC\n+\n";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_owned();

        assert!(result.is_err()); // Expect error: Incomplete record (missing quality)
        Ok(())
    }

    #[test]
    fn test_unmatched_seq_qual_length() -> Result<()> {
        let fastq_data = "@seq1 description\nATGC\n+\n!!\n";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_owned();

        assert!(result.is_err()); // Expect error: Unequal sequence and quality lengths
        Ok(())
    }

    #[test]
    fn test_edge_case_empty_data() -> Result<()> {
        let fastq_data = "";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_owned();

        assert!(result.is_ok()); // Expect error: EOF
        assert!(result.unwrap().is_none());
        Ok(())
    }

    #[test]
    fn eof_keeps_and_reprocesses_the_saved_tail_on_each_read() {
        let input = b" \t\r";
        for capacity in 1..=input.len() + 1 {
            let mut reader = FastqReader::with_capacity(capacity, input.as_slice());
            for lines in 1..=2 {
                assert!(reader.read_view().unwrap().is_none());
                assert_eq!(reader.offset(), lines);
                assert!(matches!(reader.state.state, RecordState::Empty));
            }
            let view = reader.reader.fill_view().unwrap();
            assert_eq!(view.leftover(), input);
            assert!(view.buffer().is_empty());
        }
    }

    #[test]
    fn incomplete_eof_keeps_the_parsed_record_state_and_saved_bytes() {
        let input = b"@r\nAC\n+\n";
        for capacity in 1..=input.len() + 1 {
            let mut reader = FastqReader::with_capacity(capacity, input.as_slice());
            assert!(matches!(
                reader.read_owned(),
                Err(Error::Parse {
                    line: 3,
                    source: ParseError::IncompleteRecord { .. },
                })
            ));
            let RecordState::Separator {
                header,
                sequence,
                separator,
            } = &reader.state.state
            else {
                panic!("EOF must preserve the parsed header, sequence, and separator");
            };
            assert_eq!(header, &(0..2));
            assert_eq!(sequence, &(3..5));
            assert_eq!(separator, &(6..7));
            let view = reader.reader.fill_view().unwrap();
            assert_eq!(view.leftover(), input);
            assert!(view.buffer().is_empty());
        }
    }

    #[test]
    fn fastq_eof_distinguishes_missing_and_empty_quality_lines() -> Result<()> {
        for newline in ["\n", "\r\n"] {
            let incomplete = format!("@read1{newline}{newline}+{newline}");
            let complete = format!("{incomplete}{newline}");
            for capacity in 1..=complete.len() + 1 {
                let mut reader = FastqReader::with_capacity(capacity, Cursor::new(&incomplete));
                let error = reader
                    .read_owned()
                    .expect_err("missing quality line must fail");
                assert!(matches!(
                    error,
                    Error::Parse {
                        line: 3,
                        source: ParseError::IncompleteRecord { .. }
                    }
                ));
                assert_eq!(reader.offset(), 3);

                let mut reader = FastqReader::with_capacity(capacity, Cursor::new(&complete));
                let record = reader.read_owned()?.expect("empty quality line is present");
                assert_eq!(record.id(), b"read1");
                assert!(record.sequence().is_empty());
                assert!(record.quality().is_empty());
                assert!(reader.read_owned()?.is_none());
                assert_eq!(reader.offset(), 4);
            }
        }
        Ok(())
    }

    #[test]
    fn fastq_fields_are_independent_of_crlf_buffer_boundaries() -> Result<()> {
        let input = b"@read1 description\r\nACGT\r\n+\r\nIIII\r\n";
        for capacity in 1..=input.len() + 1 {
            let mut reader = FastqReader::with_capacity(capacity, Cursor::new(input));
            let record = reader.read_owned()?.expect("missing FASTQ record");
            assert_eq!(record.id(), b"read1", "capacity {capacity}");
            assert_eq!(record.description(), Some(b"description".as_slice()));
            assert_eq!(record.sequence(), b"ACGT");
            assert_eq!(record.separator(), b"");
            assert_eq!(record.quality(), b"IIII");
            assert_eq!(reader.offset(), 4);
            assert!(reader.read_owned()?.is_none());
            assert_eq!(reader.offset(), 4);
        }
        Ok(())
    }

    #[test]
    fn skips_whitespace_only_lines_and_preserves_field_whitespace() {
        let prefix = b" \t\r\n\x0c \t\r\n";
        let first = b"@first note \t\r\n \t\r\n+ first \t\r\n \x0c\r\n";
        let between = b"\x0c\t \n";
        let second = b"@second\nAC\n+\n12\n";
        let suffix = b" \t\r\n\x0c \r";
        let input = [prefix.as_slice(), first, between, second, suffix].concat();
        for capacity in 0..=input.len() + 1 {
            let mut reader = FastqReader::with_capacity(capacity, input.as_slice());
            let record = reader.read_owned().unwrap().unwrap();
            assert_eq!(record.id(), b"first");
            assert_eq!(record.description(), Some(b"note \t".as_slice()));
            assert_eq!(record.sequence(), b" \t");
            assert_eq!(record.separator(), b" first \t");
            assert_eq!(record.quality(), b" \x0c");
            assert_eq!(record.into_parts().0, first.as_slice());
            assert_eq!(reader.offset(), 6);
            let record = reader.read_owned().unwrap().unwrap();
            assert_eq!(record.into_parts().0, second.as_slice());
            assert_eq!(reader.offset(), 11);
            for lines in 13..=14 {
                assert!(reader.read_owned().unwrap().is_none());
                assert_eq!(reader.offset(), lines);
            }
        }
    }

    #[test]
    fn a_spanning_record_does_not_copy_the_following_buffered_record() {
        let first = b"@long\r\nACGTACGTACGTACGT\r\n+\r\n1234567890123456\r\n";
        let second = b"@r\nAC\n+\n12\n";
        let input = [first.as_slice(), second].concat();
        let mut reader = FastqReader::with_capacity(first.len() - 1, input.as_slice());
        let record = reader.read_view().unwrap().unwrap();
        assert_eq!(record.into_parts().0, first);
        let buffered_pointer = reader.reader.get_mut().buffer().as_ptr();
        let record = reader.read_view().unwrap().unwrap();
        let (raw, _) = record.into_parts();
        assert_eq!(raw, second);
        assert_eq!(raw.as_ptr(), buffered_pointer);
        assert!(reader.read_view().unwrap().is_none());
    }

    #[test]
    fn owned_record_outlives_reader_after_a_borrowed_read() {
        let first = b"@first\nACGT\n+first\n1234\n";
        let second = b"@r\nTG\n+\n56\n";
        let input = [first.as_slice(), second].concat();
        let mut reader = FastqReader::with_capacity(3, input.as_slice());
        let first_record = reader.read_view().unwrap().unwrap();
        assert_eq!(first_record.into_parts().0, first);
        let record = reader.read_owned().unwrap().unwrap();
        let (raw, _) = record.into_parts();
        assert_eq!(raw, second);
        assert!(reader.read_owned().unwrap().is_none());
        drop(reader);
        assert_eq!(raw, second);
    }
}
