use std::cell::Cell;
use std::io::{self, Read, Write};
use std::rc::Rc;

use bytes::Bytes;

use biorill_fastx::{Error, FastqReader, FastqRecord, ParseError};

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
        let updated = FastqRecord::<Bytes>::try_from(parts).unwrap();
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
    record.encode_into_buf(&mut prefix);
    assert_eq!(&prefix[7..], input);
    let mut small = [0; 2];
    assert_eq!(
        record.write(&mut small.as_mut_slice()).unwrap_err().kind(),
        io::ErrorKind::WriteZero
    );
}
