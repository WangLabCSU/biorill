use std::io::{self, Write};

use biorill_fastx::{FastqField, FastqParts, FastqReader, FastqRecord};
use bytes::{BufMut, BytesMut};

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
            let encode_into = |buf: &mut dyn BufMut| {
                if use_record {
                    record.encode_into_buf(buf)
                } else {
                    fields.encode_into_buf(buf)
                }
            };
            let mut vector = b"prefix".to_vec();
            encode_into(&mut vector);
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

    let fields = FastqField::try_from(FastqParts::new("r", Some("note"), "AC", "r", "12")).unwrap();
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
