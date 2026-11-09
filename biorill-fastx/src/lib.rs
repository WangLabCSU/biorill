//! Read, edit, and write four-line FASTQ records.
//!
//! [`FastqReader`] accepts LF and CRLF, skips blank lines, and checks the record
//! format and sequence/quality lengths. Each sequence and quality must fit on
//! one line. Use `check_sequence` and `check_quality` on a record or its fields
//! to also check nucleotide and quality characters.
//!
//! Use [`FastqReader::read_view`] to borrow a record, or
//! [`FastqReader::read_owned`] to copy it into an independent record.
//! Borrowing avoids copying unless a record spans buffer refills.
//! Both methods retain the original line endings and header whitespace.
//!
//! [`FastqRecord::fields`] borrows the fields for editing. [`FastqField`] can
//! also be constructed from [`FastqParts`]. Both record types support storage
//! implementing [`stringbytes::StringBytes`].
//!
//! [`FastqRecord::write`] and [`FastqField::write`] write FASTQ with LF endings
//! and one space before a description. Separator values exclude the leading
//! `+`; an empty value produces a bare `+` line. Use `encode_into_vec` to append
//! FASTQ to a vector. The optional `bytes` feature adds `encode_into_buf` and
//! `Bytes`/`BytesMut` storage conversions.
//!
//! ```
//! use biorill_fastx::FastqReader;
//!
//! let input = b"@read1\tUMI:AC\r\nAC\r\n+\r\nII\r\n";
//! let mut reader = FastqReader::new(input.as_slice());
//! let record = reader.read_owned()?.unwrap();
//! assert_eq!(record.id(), b"read1");
//! assert_eq!(record.sequence(), b"AC");
//! assert_eq!(record.separator(), b"");
//! let mut encoded = Vec::new();
//! record.encode_into_vec(&mut encoded);
//! assert_eq!(encoded, b"@read1 UMI:AC\nAC\n+\nII\n");
//! assert_eq!(record.into_parts().0, input);
//! assert!(reader.read_owned()?.is_none());
//! # Ok::<(), biorill_fastx::Error>(())
//! ```
#![deny(unnameable_types, unreachable_pub)]

mod domain;
mod error;
mod reader;

pub use self::domain::{
    FastqField, FastqFieldError, FastqPairError, FastqParts, FastqRecord, FastqRecordError,
    LengthMismatch, ParseError,
};
pub use self::error::Error;
pub use self::reader::FastqReader;

#[cfg(test)]
mod tests {
    use std::io::{self, BufWriter, Write};

    #[cfg(feature = "bytes")]
    use bytes::{BufMut, BytesMut};

    use super::{FastqField, FastqParts, FastqReader, FastqRecord};

    pub(crate) fn valid_parts() -> FastqParts<&'static [u8]> {
        FastqParts {
            id: b"read1",
            description: Some(b"description"),
            sequence: b"AC",
            separator: b"read1",
            quality: b"12",
        }
    }

    #[test]
    fn fields_and_records_round_trip_through_a_buffered_writer() {
        let expected = b"@r note\nA\xff\n++r\n12\n";
        for capacity in [0, 3, 1024] {
            let field = FastqField::try_from(FastqParts::new(
                b"r".as_slice(),
                Some(b"note".as_slice()),
                b"A\xff",
                b"+r",
                b"12",
            ))
            .unwrap();
            let record = FastqRecord::<Vec<u8>>::from(field.clone());
            let mut writer = BufWriter::with_capacity(capacity, Vec::new());
            field.write(&mut writer).unwrap();
            record.write(&mut writer).unwrap();
            field.write(&mut writer).unwrap();
            record.write(&mut writer).unwrap();
            writer.flush().unwrap();
            assert_eq!(writer.get_ref(), &expected.repeat(4));

            let mut reader = FastqReader::new(writer.get_ref().as_slice());
            for _ in 0..4 {
                let decoded = reader.read_view().unwrap().unwrap();
                assert_eq!(decoded.sequence(), b"A\xff");
                assert_eq!(decoded.separator(), b"+r");
            }
            assert!(reader.read_view().unwrap().is_none());
        }
    }

    #[test]
    fn writing_propagates_partial_writes_and_leaves_flushing_to_the_caller() {
        let field = FastqField::try_from(FastqParts::new("r", None, "AC", "", "12")).unwrap();
        let expected = b"@r\nAC\n+\n12\n";
        for length in 0..expected.len() {
            let mut storage = vec![0; length];
            let mut slice = storage.as_mut_slice();
            let writer: &mut dyn Write = &mut slice;
            let error = field.write(writer).unwrap_err();
            assert_eq!(error.kind(), io::ErrorKind::WriteZero);
            assert_eq!(&storage, &expected[..length]);
        }

        #[derive(Default)]
        struct FailsFlush(Vec<u8>);

        impl Write for FailsFlush {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                self.0.extend_from_slice(bytes);
                Ok(bytes.len())
            }

            fn flush(&mut self) -> io::Result<()> {
                Err(io::ErrorKind::Other.into())
            }
        }

        let mut writer = BufWriter::new(FailsFlush::default());
        field.write(&mut writer).unwrap();
        assert!(writer.get_ref().0.is_empty());
        assert_eq!(writer.flush().unwrap_err().kind(), io::ErrorKind::Other);
        assert_eq!(&writer.get_ref().0, expected);
    }

    #[test]
    fn writes_records_with_and_without_descriptions() {
        for input in [
            b"@SEQ_ID desc\nACGT\n+\n1234\n".as_slice(),
            b"@SEQ_ID\nACGT\n+\n1234\n",
        ] {
            let record = FastqReader::new(input).read_owned().unwrap().unwrap();
            let mut output = Vec::new();
            record.write(&mut output).unwrap();
            assert_eq!(output, input);
        }
    }

    #[test]
    fn fields_and_records_encode_raw_bytes_into_writers_and_buf_mut_targets() {
        for (description, expected) in [
            (None, b"@read1\nA\xff\n++read1\n12\n".as_slice()),
            (Some(b"".as_slice()), b"@read1 \nA\xff\n++read1\n12\n"),
            (
                Some(b"note".as_slice()),
                b"@read1 note\nA\xff\n++read1\n12\n",
            ),
        ] {
            let fields = FastqField::try_from(FastqParts::new(
                b"read1".as_slice(),
                description,
                b"A\xff",
                b"+read1",
                b"12",
            ))
            .unwrap();
            let record = FastqRecord::<Vec<u8>>::from(fields.clone());
            assert_eq!(
                {
                    let mut encoded = Vec::new();
                    fields.encode_into_vec(&mut encoded);
                    encoded
                },
                expected
            );
            assert_eq!(
                {
                    let mut encoded = Vec::new();
                    record.encode_into_vec(&mut encoded);
                    encoded
                },
                expected
            );

            for use_record in [false, true] {
                let mut vector = b"prefix".to_vec();
                if use_record {
                    record.encode_into_vec(&mut vector);
                } else {
                    fields.encode_into_vec(&mut vector);
                }
                assert_eq!(vector, [b"prefix".as_slice(), expected].concat());

                let mut reserved = Vec::with_capacity(vector.len() * 2);
                reserved.extend_from_slice(b"prefix");
                let pointer = reserved.as_ptr();
                if use_record {
                    record.encode_into_vec(&mut reserved);
                    record.encode_into_vec(&mut reserved);
                } else {
                    fields.encode_into_vec(&mut reserved);
                    fields.encode_into_vec(&mut reserved);
                }
                assert_eq!(
                    reserved,
                    [b"prefix".as_slice(), expected, expected].concat()
                );
                assert_eq!(reserved.as_ptr(), pointer);

                #[cfg(feature = "bytes")]
                {
                    let encode_into = |buf: &mut dyn BufMut| {
                        if use_record {
                            record.encode_into_buf(buf)
                        } else {
                            fields.encode_into_buf(buf)
                        }
                    };
                    let mut appended = b"prefix".to_vec();
                    encode_into(&mut appended);
                    assert_eq!(appended, vector);

                    let mut bytes = BytesMut::from(b"prefix".as_slice());
                    encode_into(&mut bytes);
                    assert_eq!(bytes.as_ref(), vector);

                    let mut fixed = vec![0; expected.len()];
                    let mut target = fixed.as_mut_slice();
                    encode_into(&mut target);
                    assert!(target.is_empty());
                    assert_eq!(fixed, expected);

                    let mut first = [0; 3];
                    let mut second = vec![0; expected.len() - first.len()];
                    encode_into(&mut first.as_mut_slice().chain_mut(second.as_mut_slice()));
                    assert_eq!([first.as_slice(), second.as_slice()].concat(), expected);
                }

                let mut written = b"prefix".to_vec();
                let writer: &mut dyn Write = &mut written;
                if use_record {
                    record.write(writer)
                } else {
                    fields.write(writer)
                }
                .unwrap();
                assert_eq!(written, vector);
            }
        }
    }

    #[test]
    fn writes_retry_interruptions_and_short_writes_and_stop_on_io_errors() {
        struct LimitedWriter {
            output: Vec<u8>,
            limit: usize,
            interrupted: bool,
            failed: bool,
        }

        impl Write for LimitedWriter {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                assert!(!self.failed, "writing continued after an error");
                if !self.interrupted {
                    self.interrupted = true;
                    return Err(io::ErrorKind::Interrupted.into());
                }
                if self.output.len() == self.limit {
                    self.failed = true;
                    return Err(io::ErrorKind::Other.into());
                }
                let count = bytes.len().min(2).min(self.limit - self.output.len());
                self.output.extend_from_slice(&bytes[..count]);
                Ok(count)
            }

            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let fields =
            FastqField::try_from(FastqParts::new("r", Some("note"), "AC", "r", "12")).unwrap();
        let record = FastqRecord::<Vec<u8>>::from(fields.clone());
        let expected = b"@r note\nAC\n+r\n12\n";
        for use_record in [false, true] {
            for limit in 0..=expected.len() {
                let mut writer = LimitedWriter {
                    output: Vec::new(),
                    limit,
                    interrupted: false,
                    failed: false,
                };
                let result = if use_record {
                    record.write(&mut writer)
                } else {
                    fields.write(&mut writer)
                };
                assert_eq!(result.is_ok(), limit == expected.len());
                assert_eq!(writer.output, expected[..limit]);
                assert_eq!(writer.failed, limit < expected.len());
            }
        }
    }

    #[cfg(feature = "bytes")]
    #[test]
    fn encoding_into_a_full_buffer_panics_after_an_optional_prefix() {
        let fields = FastqField::try_from(FastqParts::new("r", None, "AC", "", "12")).unwrap();
        let record = FastqRecord::<Vec<u8>>::from(fields.clone());
        let expected = b"@r\nAC\n+\n12\n";
        for use_record in [false, true] {
            for capacity in 0..expected.len() {
                let mut storage = vec![0; capacity];
                let mut target = storage.as_mut_slice();
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    if use_record {
                        record.encode_into_buf(&mut target)
                    } else {
                        fields.encode_into_buf(&mut target)
                    }
                }));
                assert!(result.is_err());
                let written = capacity - target.len();
                assert_eq!(storage[..written], expected[..written]);
            }
        }
    }

    #[test]
    fn records_survive_later_reads_and_round_trip_across_all_buffer_sizes() {
        let input = b" \t\r\n@r1\tMIRE{UMI:AC}\r\nA\xff\r\n+r1\r\n!~\r\n\n@r2\nTG\n+\n12";
        for capacity in 0..=input.len() + 1 {
            let mut reader = FastqReader::with_capacity(capacity, input.as_slice());
            let first = reader.read_owned().unwrap().unwrap();
            let second = reader.read_owned().unwrap().unwrap();
            for _ in 0..3 {
                assert!(reader.read_owned().unwrap().is_none());
                assert_eq!(reader.offset(), 10);
            }
            assert_eq!(first.id(), b"r1");
            assert_eq!(first.description().unwrap(), b"MIRE{UMI:AC}");
            assert_eq!(first.sequence(), b"A\xff");
            assert_eq!(first.separator(), b"r1");
            assert_eq!(first.quality(), b"!~");
            assert_eq!(
                {
                    let mut encoded = Vec::new();
                    FastqRecord::<&[u8]>::from(&first).encode_into_vec(&mut encoded);
                    encoded
                },
                b"@r1 MIRE{UMI:AC}\nA\xff\n+r1\n!~\n"
            );
            let expected = b"@r2\nTG\n+\n12\n";
            assert_eq!(
                {
                    let mut encoded = Vec::new();
                    second.encode_into_vec(&mut encoded);
                    encoded
                },
                expected
            );
            let mut parts = first.parts();
            parts.sequence = second.sequence();
            parts.quality = second.quality();
            parts.description = Some(b"changed");
            let updated = FastqRecord::<Vec<u8>>::try_from(parts).unwrap();
            let (data, ranges) = second.into_parts();
            assert_eq!(data.as_slice(), b"@r2\nTG\n+\n12");
            assert_eq!(&data[ranges.sequence], b"TG");
            assert_eq!(
                {
                    let mut encoded = Vec::new();
                    updated.encode_into_vec(&mut encoded);
                    encoded
                },
                b"@r1 changed\nTG\n+r1\n12\n"
            );
            assert_eq!(
                FastqRecord::<&[u8]>::from(&first).into_parts().0,
                b"@r1\tMIRE{UMI:AC}\r\nA\xff\r\n+r1\r\n!~\r\n"
            );
        }
    }

    #[test]
    fn separator_fields_skip_one_marker_and_round_trip_all_remaining_bytes() {
        for content in [
            b"".as_slice(),
            b"read1",
            b"+read1",
            b"++",
            b" \tread1 \t",
            b"\xff",
        ] {
            let expected = [b"@read1\nAC\n+".as_slice(), content, b"\n12\n"].concat();
            for newline in [b"\n".as_slice(), b"\r\n"] {
                let input = [
                    b"@read1".as_slice(),
                    newline,
                    b"AC",
                    newline,
                    b"+",
                    content,
                    newline,
                    b"12",
                    newline,
                ]
                .concat();
                for capacity in 0..=input.len() + 1 {
                    let mut reader = FastqReader::with_capacity(capacity, input.as_slice());
                    let record = reader.read_owned().unwrap().unwrap();
                    assert_eq!(record.separator(), content);
                    assert_eq!(record.parts().separator, content);
                    assert_eq!(*record.fields().separator(), content);
                    assert_eq!(
                        FastqRecord::<&[u8]>::from(&record).into_parts().0,
                        input.as_slice()
                    );
                    assert_eq!(
                        {
                            let mut encoded = Vec::new();
                            record.encode_into_vec(&mut encoded);
                            encoded
                        },
                        expected
                    );
                    let mut written = Vec::new();
                    record.write(&mut written).unwrap();
                    assert_eq!(written, expected);
                    let rebuilt = FastqRecord::<Vec<u8>>::try_from(record.parts()).unwrap();
                    assert_eq!(
                        FastqRecord::<&[u8]>::from(&rebuilt).into_parts().0,
                        [b"read1AC".as_slice(), content, b"12"].concat()
                    );
                    assert_eq!(
                        {
                            let mut encoded = Vec::new();
                            rebuilt.encode_into_vec(&mut encoded);
                            encoded
                        },
                        expected
                    );
                    assert_eq!(rebuilt.separator(), content);
                    let (data, ranges) = record.into_parts();
                    let separator_start =
                        b"@read1".len() + newline.len() + b"AC".len() + newline.len() + 1;
                    assert_eq!(
                        ranges.separator,
                        separator_start..separator_start + content.len()
                    );
                    assert_eq!(data[ranges.separator.start - 1], b'+');
                    assert_eq!(&data[ranges.separator], content);
                    assert!(reader.read_owned().unwrap().is_none());
                }
            }
        }
    }

    #[test]
    fn writing_handles_short_writes_and_propagates_output_errors() {
        struct ShortWriter(Vec<u8>);
        impl Write for ShortWriter {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                let count = bytes.len().min(2);
                self.0.extend_from_slice(&bytes[..count]);
                Ok(count)
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let input = b"@read1\nAC\n+read1\n12\n";
        let record = FastqReader::new(input.as_slice())
            .read_owned()
            .unwrap()
            .unwrap();
        let mut writer = ShortWriter(Vec::new());
        record.write(&mut writer).unwrap();
        assert_eq!(writer.0, input);
        let mut prefix = b"prefix\n".to_vec();
        record.encode_into_vec(&mut prefix);
        assert_eq!(&prefix[7..], input);
        let mut small = [0; 2];
        assert_eq!(
            record.write(&mut small.as_mut_slice()).unwrap_err().kind(),
            io::ErrorKind::WriteZero
        );
    }

    #[test]
    fn text_display_matches_fastq_writer_for_optional_descriptions_and_empty_fields() {
        for description in [None, Some(""), Some("说明\ttext")] {
            for (id, sequence, separator, quality) in [
                ("读1", "AC", "读1", "12"),
                ("", "", "", ""),
                ("read1", "é", "read1", "12"),
                ("read1", "AC", "+read1", "12"),
            ] {
                let borrowed = FastqField::try_from(FastqParts {
                    id,
                    description,
                    sequence,
                    separator,
                    quality,
                })
                .unwrap();
                let owned = FastqField::<String>::try_from_fields(borrowed.clone()).unwrap();
                let record = FastqRecord::<Vec<u8>>::from(borrowed.clone());
                let mut expected = Vec::new();
                record.write(&mut expected).unwrap();

                let record_fields = record.fields();
                let record_text = record_fields.try_str_view().unwrap();
                assert_eq!(record_text.description().is_some(), description.is_some());
                assert_eq!(record_text.to_string().as_bytes(), expected);
                assert_eq!(borrowed.to_string().as_bytes(), expected);
                assert_eq!(owned.to_string().as_bytes(), expected);
                assert_eq!(owned.as_deref().to_string().as_bytes(), expected);

                let bytes = FastqField::<Vec<u8>>::from_fields(owned);
                let checked = bytes.try_str_view().unwrap();
                assert_eq!(checked.to_string().as_bytes(), expected);
                assert_eq!(checked.sequence().as_ptr(), bytes.sequence().as_ptr());
            }
        }
    }
}
