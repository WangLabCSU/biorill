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
    use std::cell::Cell;
    use std::collections::VecDeque;
    use std::io::{self, Cursor, Read};
    use std::rc::Rc;

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

    #[test]
    fn views_borrow_records_directly_from_the_read_destination() {
        struct TrackedInput<'a> {
            data: &'a [u8],
            destination: Rc<Cell<usize>>,
            reads: Rc<Cell<usize>>,
        }

        impl Read for TrackedInput<'_> {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.reads.get() == 0 {
                    self.destination.set(buf.as_ptr() as usize);
                }
                self.reads.set(self.reads.get() + 1);
                self.data.read(buf)
            }
        }

        let prefix = b" \t\r\n";
        let first = b"@first\tnote\r\nAC\r\n+first\r\n12\r\n";
        let second = b"@second\nTG\n+\n34\n";
        let input = [prefix.as_slice(), first, second].concat();
        let destination = Rc::new(Cell::new(0));
        let reads = Rc::new(Cell::new(0));
        let mut reader = FastqReader::with_capacity(
            input.len() + 16,
            TrackedInput {
                data: &input,
                destination: Rc::clone(&destination),
                reads: Rc::clone(&reads),
            },
        );
        let record1 = reader.read_view().unwrap().unwrap();
        let (raw1, _) = record1.into_parts();
        assert_eq!(raw1.as_ptr() as usize, destination.get() + prefix.len());
        assert_eq!(raw1, first.as_slice());
        let record2 = reader.read_view().unwrap().unwrap();
        let (raw2, _) = record2.into_parts();
        assert_eq!(
            raw2.as_ptr() as usize,
            destination.get() + prefix.len() + first.len()
        );
        assert_eq!(raw2, second.as_slice());
        assert_eq!(reads.get(), 1);
        for _ in 0..3 {
            assert!(reader.read_view().unwrap().is_none());
        }
    }

    #[test]
    fn large_records_span_many_reads_without_changing_retained_records() {
        let first = b"@first\nAC\n+\n12\n";
        let last = b"@last\nTG\n+\n34\n";
        let sequence = vec![b'A'; 8192];
        let quality = vec![b'!'; sequence.len()];
        let middle = [
            b"@large\r\n".as_slice(),
            &sequence,
            b"\r\n+large\r\n",
            &quality,
            b"\r\n",
        ]
        .concat();
        let input = [first.as_slice(), &middle, last].concat();
        for capacity in [1, 7, 128] {
            let mut reader = FastqReader::with_capacity(capacity, input.as_slice());
            let record1 = reader.read_owned().unwrap().unwrap();
            let record2 = reader.read_owned().unwrap().unwrap();
            let record3 = reader.read_owned().unwrap().unwrap();
            assert!(reader.read_owned().unwrap().is_none());
            assert_eq!(record2.sequence(), sequence);
            assert_eq!(record2.quality(), quality);
            assert_eq!(record1.into_parts().0, first.as_slice());
            assert_eq!(record2.into_parts().0, middle);
            assert_eq!(record3.into_parts().0, last.as_slice());
        }
    }

    #[test]
    fn final_cr_is_content_and_line_markers_in_quality_are_not_record_boundaries() {
        for (input, sequence, quality) in [
            (
                b"@read\nAC\n+\n1\r".as_slice(),
                b"AC".as_slice(),
                b"1\r".as_slice(),
            ),
            (b"@read\r\nAC\r\n+\r\n@+\r\n", b"AC", b"@+"),
        ] {
            for capacity in 1..=input.len() {
                let mut reader = FastqReader::with_capacity(capacity, input);
                let record = reader.read_owned().unwrap().unwrap();
                assert_eq!(record.sequence(), sequence);
                assert_eq!(record.quality(), quality);
                assert_eq!(record.into_parts().0, input);
                assert!(reader.read_owned().unwrap().is_none());
            }
        }
    }

    #[test]
    fn carriage_returns_in_a_saved_prefix_do_not_trim_unrelated_bytes() {
        for (input, id) in [
            (b"@r\rx\nAC\n+\n12\n".as_slice(), b"r\rx".as_slice()),
            (b"@long_identifier\r\nAC\n+\n12\n", b"long_identifier"),
        ] {
            for capacity in 0..=input.len() + 1 {
                let mut reader = FastqReader::with_capacity(capacity, input);
                let record = reader.read_view().unwrap().unwrap();
                assert_eq!(record.id(), id, "capacity {capacity}");
                assert_eq!(record.sequence(), b"AC");
                assert_eq!(record.quality(), b"12");
                assert_eq!(record.into_parts().0, input);
                assert_eq!(reader.offset(), 4);
                assert!(reader.read_view().unwrap().is_none());
            }
        }
    }

    #[test]
    fn invalid_header_and_separator_are_reported_before_a_later_read_failure() {
        struct FailsNext<'a>(&'a [u8]);

        impl Read for FailsNext<'_> {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.0.is_empty() {
                    Err(io::ErrorKind::Other.into())
                } else {
                    self.0.read(buf)
                }
            }
        }

        for (input, expected_line) in [(b"read\n".as_slice(), 1), (b"@read\nAC\nwrong\n", 3)] {
            for capacity in 0..=input.len() + 1 {
                for borrowed in [false, true] {
                    let mut reader = FastqReader::with_capacity(capacity, FailsNext(input));
                    let error = if borrowed {
                        reader.read_view().unwrap_err()
                    } else {
                        reader.read_owned().unwrap_err()
                    };
                    let Error::Parse { line, source } = error else {
                        panic!("invalid input must be rejected before requesting more bytes");
                    };
                    assert_eq!(line, expected_line);
                    assert!(matches!(
                        source,
                        ParseError::InvalidHead { .. } | ParseError::InvalidSep { .. }
                    ));
                }
            }
        }
    }

    #[test]
    fn final_record_without_newline_does_not_read_again_after_eof() {
        struct FailsAfterEof<'a> {
            data: &'a [u8],
            eof: bool,
        }

        impl Read for FailsAfterEof<'_> {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if !self.data.is_empty() {
                    self.data.read(buf)
                } else if self.eof {
                    Err(io::ErrorKind::Other.into())
                } else {
                    self.eof = true;
                    Ok(0)
                }
            }
        }

        let input = b"@r note\r\nAC\r\n+r\r\n12";
        for capacity in 0..=input.len() + 1 {
            for borrowed in [false, true] {
                let mut reader = FastqReader::with_capacity(
                    capacity,
                    FailsAfterEof {
                        data: input,
                        eof: false,
                    },
                );
                if borrowed {
                    let record = reader.read_view().unwrap().unwrap();
                    assert_eq!(record.into_parts().0, input);
                } else {
                    let record = reader.read_owned().unwrap().unwrap();
                    assert_eq!(record.into_parts().0, input);
                }
                assert_eq!(reader.offset(), 4);
                let Error::Read { source } = reader.read_owned().unwrap_err() else {
                    panic!("the read after EOF must be left to the next call");
                };
                assert_eq!(source.kind(), io::ErrorKind::Other);
            }
        }
    }

    #[test]
    fn error_messages_preserve_utf8_across_chunk_boundaries() {
        for (input, expected, expected_line) in [
            ("réad\n", "expected '@' at record start\nréad", 1),
            (
                "@réad 描述\nAC\n错误\n",
                "expected '+' at separator line start\n@réad 描述\nAC\n错误",
                3,
            ),
            (
                "@réad 描述\nAC\n+réad\n1\n",
                "sequence and quality lengths do not match (2 vs 1)\n@réad 描述\nAC\n+réad\n1",
                4,
            ),
            ("@réad 描述", "incomplete record\nréad描述", 1),
        ] {
            for capacity in 0..=input.len() + 1 {
                let mut reader = FastqReader::with_capacity(capacity, input.as_bytes());
                let Error::Parse { line, source } = reader.read_view().unwrap_err() else {
                    panic!("expected a parse error");
                };
                assert_eq!(line, expected_line);
                assert_eq!(source.to_string(), expected, "capacity {capacity}");
            }
        }
    }

    #[test]
    fn failures_identify_the_kind_and_physical_line_across_buffer_sizes() {
        let cases: &[(&[u8], &str, usize)] = &[
            (b"\nread1\nAC\n+\n12\n", "header", 2),
            (b"@read1\n", "incomplete", 1),
            (b"@read1\nAC\n", "incomplete", 2),
            (b"@read1\nAC\n+\n", "incomplete", 3),
            (b"@read1\nAC\nwrong\n12\n", "separator", 3),
            (b"@read1\nAC\n\n12\n", "separator", 3),
            (b"\r\n@read1\r\nAC\r\n+\r\n1\r\n", "length", 5),
        ];
        for &(input, expected_kind, expected_line) in cases {
            for capacity in 0..=input.len() + 1 {
                let mut reader = FastqReader::with_capacity(capacity, input);
                let error = reader.read_owned().unwrap_err();
                let Error::Parse {
                    line,
                    source: parse,
                } = error
                else {
                    panic!("expected a parse error: {error}");
                };
                let kind = match parse {
                    ParseError::InvalidHead { .. } => "header",
                    ParseError::InvalidSep { .. } => "separator",
                    ParseError::IncompleteRecord { .. } => "incomplete",
                    ParseError::InvalidRecord { .. } => "record",
                    ParseError::UnequalLength { seq, qual, .. } => {
                        assert_eq!((seq, qual), (2, 1));
                        "length"
                    }
                };
                assert_eq!((kind, line), (expected_kind, expected_line));
                assert_eq!(reader.offset(), expected_line);
            }
        }
    }

    #[test]
    fn input_failures_inside_each_record_line_remain_io_errors() {
        struct FailsAfter<'a>(&'a [u8]);
        impl Read for FailsAfter<'_> {
            fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
                if self.0.is_empty() {
                    Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "truncated input",
                    ))
                } else {
                    self.0.read(out)
                }
            }
        }
        let input = b"@read1 description\nAC\n+read1\n12\n";
        for split in 0..input.len() {
            let mut reader = FastqReader::with_capacity(3, FailsAfter(&input[..split]));
            let error = reader.read_owned().unwrap_err();
            let Error::Read { source } = error else {
                panic!("I/O failures must not become EOF or parse failures: {error}");
            };
            assert_eq!(
                reader.offset() + 1,
                input[..split].iter().filter(|&&byte| byte == b'\n').count() + 1
            );
            assert_eq!(source.kind(), io::ErrorKind::UnexpectedEof);
            assert_eq!(source.to_string(), "truncated input");
        }
    }

    #[test]
    fn reading_modes_can_alternate_and_owned_records_survive_later_reads() {
        let first = b"@first\tnote\r\nA\xff\r\n++first\r\n12\r\n";
        let second = b"@second\nTG\n+\n34\n";
        let third = b"@third\n\n+third\n\n";
        let last = b"@last\nAC\n+\n56";
        let input = [b" \t\n".as_slice(), first, second, third, last].concat();
        for capacity in 0..=input.len() + 1 {
            let mut reader = FastqReader::with_capacity(capacity, input.as_slice());
            let view = reader.read_view().unwrap().unwrap();
            assert_eq!(view.id(), b"first");
            assert_eq!(view.sequence(), b"A\xff");
            assert_eq!(view.into_parts().0, first);

            let owned = reader.read_owned().unwrap().unwrap();
            let vector = reader.read_owned().unwrap().unwrap();
            assert_eq!(reader.read_view().unwrap().unwrap().into_parts().0, last);
            assert!(reader.read_owned().unwrap().is_none());
            assert_eq!(reader.offset(), 17);
            drop(reader);
            assert_eq!(owned.into_parts().0, second.as_slice());
            assert_eq!(vector.into_parts().0, third.as_slice());
        }
    }

    #[test]
    fn input_can_resume_between_records_after_eof_with_a_partial_header() {
        struct Chunks<'a>(VecDeque<Option<&'a [u8]>>);

        impl Read for Chunks<'_> {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                loop {
                    match self.0.front_mut() {
                        None => return Ok(0),
                        Some(None) => {
                            self.0.pop_front();
                            return Ok(0);
                        }
                        Some(Some(bytes)) if bytes.is_empty() => {
                            self.0.pop_front();
                        }
                        Some(Some(bytes)) => return bytes.read(buf),
                    }
                }
            }
        }

        for capacity in [0, 1, 7, 64] {
            let source = Chunks(VecDeque::from([
                Some(b"@first\nA\n+\n1".as_slice()),
                None,
                None,
                Some(b"@sec"),
                Some(b"ond\nTG\n+\n23\n"),
            ]));
            let mut reader = FastqReader::with_capacity(capacity, source);
            let first = reader.read_owned().unwrap().unwrap();
            assert!(reader.read_owned().unwrap().is_none());
            assert_eq!(reader.offset(), 4);
            let second = reader.read_view().unwrap().unwrap();
            assert_eq!(second.id(), b"second");
            assert_eq!(second.sequence(), b"TG");
            assert_eq!(second.quality(), b"23");
            assert_eq!(reader.offset(), 8);
            assert!(reader.read_owned().unwrap().is_none());
            assert_eq!(first.into_parts().0, b"@first\nA\n+\n1".as_slice());
        }
    }
}
