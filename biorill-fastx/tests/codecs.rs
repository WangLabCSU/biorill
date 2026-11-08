use std::collections::VecDeque;
use std::io::{self, BufWriter, Read, Write};

use biorill_fastx::{FastqField, FastqParts, FastqReader, FastqRecord};

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
