use std::cell::Cell;
use std::ops::Range;
use std::rc::Rc;
use std::sync::Arc;

use biorill_fastx::{FastqField, FastqParts, FastqRecord, FastqRecordError};
use bytes::{Bytes, BytesMut};
use stringbytes::StringBytes;

const DATA: &[u8] = b"12|AC|note|read1|+|tail";
const FASTQ: &[u8] = b"@read1 note\nAC\n+\n12\n";

#[derive(Clone, Debug, PartialEq, Eq)]
struct Custom(Vec<u8>);

impl AsRef<[u8]> for Custom {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

fn ranges() -> FastqParts<Range<usize>> {
    FastqParts {
        id: 11..16,
        description: Some(6..10),
        sequence: 3..5,
        separator: 18..18,
        quality: 0..2,
    }
}

fn check_storage<T: StringBytes + AsRef<[u8]>>(data: T) {
    let pointer = data.as_ref().as_ptr();
    let record = FastqRecord::<T>::try_from((data, ranges())).unwrap();
    assert_eq!(
        FastqRecord::<&[u8]>::from(&record).into_parts().0.as_ptr(),
        pointer
    );
    check_record(record);
}

fn unpack_record<T>(record: FastqRecord<T>) -> (T, FastqParts<Range<usize>>) {
    record.into_parts()
}

fn unpack_fields<T>(fields: FastqField<T>) -> FastqParts<T> {
    fields.into_parts()
}

fn check_fields<T: StringBytes + AsRef<[u8]>>(fields: FastqField<T>) {
    assert_eq!(fields.id().as_ref(), b"read1");
    assert_eq!(
        fields.description().map(AsRef::as_ref),
        Some(b"note".as_slice())
    );
    assert_eq!(fields.sequence().as_ref(), b"AC");
    assert_eq!(fields.separator().as_ref(), b"");
    assert_eq!(fields.quality().as_ref(), b"12");
    let borrowed = fields.as_ref();
    assert_eq!(unpack_fields(borrowed).sequence.as_ref(), b"AC");
    assert_eq!(
        {
            let mut encoded = Vec::new();
            FastqRecord::<Bytes>::from(fields.as_ref()).encode_into_vec(&mut encoded);
            encoded
        },
        FASTQ
    );
    assert_eq!(
        {
            let mut encoded = Vec::new();
            FastqRecord::<BytesMut>::from(fields.as_ref()).encode_into_vec(&mut encoded);
            encoded
        },
        FASTQ
    );
    assert_eq!(
        {
            let mut encoded = Vec::new();
            FastqRecord::<Box<[u8]>>::from(fields.as_ref()).encode_into_vec(&mut encoded);
            encoded
        },
        FASTQ
    );
    assert_eq!(
        {
            let mut encoded = Vec::new();
            FastqRecord::<Vec<u8>>::from(fields).encode_into_vec(&mut encoded);
            encoded
        },
        FASTQ
    );
}

fn check_record<T: StringBytes + AsRef<[u8]>>(record: FastqRecord<T>) {
    let pointer = FastqRecord::<&[u8]>::from(&record).into_parts().0.as_ptr();
    assert_eq!(FastqRecord::<&[u8]>::from(&record).into_parts().0, DATA);
    assert_eq!(record.id(), b"read1");
    assert_eq!(record.description(), Some(b"note".as_slice()));
    assert_eq!(record.sequence(), b"AC");
    assert_eq!(record.separator(), b"");
    assert_eq!(record.quality(), b"12");
    assert_eq!(record.parts().sequence, b"AC");
    check_fields(record.fields());
    let borrowed: FastqRecord<&[u8]> = (&record).into();
    let (borrowed, borrowed_ranges) = borrowed.into_parts();
    assert_eq!(borrowed.as_ptr(), pointer);
    assert_eq!(borrowed_ranges, ranges());
    assert_eq!(
        {
            let mut encoded = Vec::new();
            record.encode_into_vec(&mut encoded);
            encoded
        },
        FASTQ
    );
    let mut written = Vec::new();
    record.write(&mut written).unwrap();
    assert_eq!(written, FASTQ);
    let mut appended = b"prefix".to_vec();
    record.encode_into_buf(&mut appended);
    assert_eq!(appended, [b"prefix".as_slice(), FASTQ].concat());
    let (data, returned_ranges) = unpack_record(record);
    assert_eq!(data.as_ref().as_ptr(), pointer);
    assert_eq!(data.as_ref(), DATA);
    assert_eq!(returned_ranges, ranges());
}

#[test]
fn standard_storage_and_borrowed_custom_bytes_retain_layout_and_encode_fields() {
    check_storage(DATA);
    check_storage(DATA.to_vec());
    check_storage(Bytes::from_static(DATA));
    check_storage(BytesMut::from(DATA));
    let custom = Custom(DATA.to_vec());
    check_storage(custom.as_ref());
}

#[test]
fn stringbytes_wrappers_borrow_record_bytes_and_preserve_ranges() {
    fn check<T: StringBytes>(data: T) {
        let pointer = data.as_bytes().as_ptr();
        let record = FastqRecord::try_from((data, ranges())).unwrap();
        assert_eq!(record.id(), b"read1");
        assert_eq!(record.sequence(), b"AC");

        let borrowed = FastqRecord::<&[u8]>::from(&record);
        let (bytes, byte_ranges) = borrowed.into_parts();
        assert_eq!(bytes.as_ptr(), pointer);
        assert_eq!(bytes, DATA);
        assert_eq!(byte_ranges, ranges());

        let text = FastqRecord::<&str>::try_from(&record).unwrap();
        assert_eq!(text.parts(), record.parts());
        assert_eq!(text.into_parts().0.as_ptr(), pointer);
        let mut encoded = Vec::new();
        record.encode_into_vec(&mut encoded);
        assert_eq!(encoded, FASTQ);

        let (data, returned_ranges) = record.into_parts();
        assert_eq!(data.as_bytes().as_ptr(), pointer);
        assert_eq!(returned_ranges, ranges());
    }

    check(Rc::new(String::from_utf8(DATA.to_vec()).unwrap()));
    check(Arc::new(Box::new(DATA.to_vec())));
    check(Rc::new(Arc::new(Box::new(DATA.to_vec()))));
    check(Box::new(DATA.to_vec()) as Box<dyn StringBytes>);
    check(Rc::new(DATA.to_vec()) as Rc<dyn StringBytes>);
    check(Arc::new(String::from_utf8(DATA.to_vec()).unwrap()) as Arc<dyn StringBytes>);
    check(DATA.to_vec().as_mut_slice());
}

#[test]
fn as_ref_borrows_the_buffer_without_calling_custom_as_ref() {
    // A local trait object is still StringBytes, but a custom AsRef may expose
    // a different byte sequence. Record borrowing must never call that impl.
    trait LocalBytes: StringBytes {}
    impl LocalBytes for Vec<u8> {}
    impl AsRef<dyn LocalBytes> for Vec<u8> {
        fn as_ref(&self) -> &(dyn LocalBytes + 'static) {
            static OTHER: Vec<u8> = Vec::new();
            &OTHER
        }
    }

    // Borrowing the actual T requires no trait bound or user-defined conversion.
    fn borrow<T>(record: &FastqRecord<T>) -> FastqRecord<&T> {
        record.as_ref()
    }

    let data = DATA.to_vec();
    assert!(AsRef::<dyn LocalBytes>::as_ref(&data).as_bytes().is_empty());
    let pointer = data.as_ptr();
    let record = FastqRecord::try_from((data, ranges())).unwrap();
    let borrowed: FastqRecord<&Vec<u8>> = borrow(&record);
    assert_eq!(borrowed.parts(), record.parts());
    let (buffer, returned_ranges) = borrowed.into_parts();
    assert_eq!(buffer.as_ptr(), pointer);
    assert_eq!(buffer, DATA);
    assert_eq!(returned_ranges, ranges());
}

#[test]
fn generic_borrows_retain_storage_and_ranges_without_revalidation() {
    let data = Rc::new(String::from("12|é|note|read1|+|tail"));
    let pointer = data.as_ptr();
    let record = FastqRecord::try_from((data, ranges())).unwrap();

    // Borrow the Rc itself, then use the explicit text conversion.
    let view: FastqRecord<&Rc<String>> = record.as_ref();
    let text = FastqRecord::<&str>::try_from(&view).unwrap();
    assert_eq!(text.sequence(), "é".as_bytes());
    assert_eq!(text.quality(), b"12");
    let (text, text_ranges) = text.into_parts();
    assert_eq!(text.as_ptr(), pointer);
    assert_eq!(text_ranges, ranges());

    let (bytes, byte_ranges) = FastqRecord::<&[u8]>::from(&view).into_parts();
    assert_eq!(bytes.as_ptr(), pointer);
    assert_eq!(bytes, text.as_bytes());
    assert_eq!(byte_ranges, text_ranges);
}

#[test]
fn fields_borrowed_from_custom_parts_support_checked_updates() {
    let parts = FastqParts::new(
        Custom(b"read0".to_vec()),
        None,
        Custom(b"A".to_vec()),
        Custom(b"read0".to_vec()),
        Custom(b"1".to_vec()),
    );
    let mut fields = FastqField::try_from(FastqParts::new(
        parts.id.as_ref(),
        parts.description.as_ref().map(AsRef::as_ref),
        parts.sequence.as_ref(),
        parts.separator.as_ref(),
        parts.quality.as_ref(),
    ))
    .unwrap();
    let original = fields.clone();
    assert!(fields.use_sequence_and_quality(b"AC", b"1").is_err());
    assert_eq!(fields, original);

    fields.use_id(b"read1").unwrap();
    fields.use_description(Some(b"note")).unwrap();
    fields.use_separator(b"").unwrap();
    fields.use_sequence_and_quality(b"AC", b"12").unwrap();
    check_fields(fields);
}

#[test]
fn owned_record_copies_borrowed_parts_from_custom_storage() {
    struct Data {
        bytes: &'static [u8],
        empty: Rc<Cell<bool>>,
    }

    impl AsRef<[u8]> for Data {
        fn as_ref(&self) -> &[u8] {
            if self.empty.get() {
                b""
            } else {
                self.bytes
            }
        }
    }

    let empty = Rc::new(Cell::new(false));
    let data = |bytes| Data {
        bytes,
        empty: Rc::clone(&empty),
    };
    let parts = FastqParts::new(
        data(b"read1"),
        Some(data(b"note")),
        data(b"AC"),
        data(b""),
        data(b"12"),
    );
    let record = FastqRecord::<Bytes>::try_from(FastqParts::new(
        parts.id.as_ref(),
        parts.description.as_ref().map(AsRef::as_ref),
        parts.sequence.as_ref(),
        parts.separator.as_ref(),
        parts.quality.as_ref(),
    ))
    .unwrap();

    empty.set(true);
    assert_eq!(
        FastqRecord::<&[u8]>::from(&record).into_parts().0,
        b"read1noteAC12"
    );
    assert_eq!(
        {
            let mut encoded = Vec::new();
            record.encode_into_vec(&mut encoded);
            encoded
        },
        FASTQ
    );
    assert_eq!(record.sequence(), b"AC");
    assert_eq!(record.quality(), b"12");
}

#[test]
fn conversion_to_bytes_retains_unformatted_bytes_and_ranges() {
    let record = FastqRecord::<BytesMut>::try_from((BytesMut::from(DATA), ranges())).unwrap();
    let pointer = FastqRecord::<&[u8]>::from(&record).into_parts().0.as_ptr();
    let frozen = FastqRecord::<Bytes>::from(record);
    assert_eq!(
        FastqRecord::<&[u8]>::from(&frozen).into_parts().0.as_ptr(),
        pointer
    );
    assert_eq!(FastqRecord::<&[u8]>::from(&frozen).into_parts().0, DATA);
    assert_eq!(
        {
            let mut encoded = Vec::new();
            frozen.encode_into_vec(&mut encoded);
            encoded
        },
        FASTQ
    );
    assert_eq!(frozen.into_parts().1, ranges());
}

#[test]
fn from_conversions_reuse_unique_storage_and_preserve_field_ranges() {
    for spare_capacity in [0, 32] {
        let mut data = Vec::with_capacity(DATA.len() + spare_capacity);
        data.extend_from_slice(DATA);
        let pointer = data.as_ptr();
        let vector = FastqRecord::<Vec<u8>>::try_from((data, ranges())).unwrap();

        let frozen: FastqRecord<Bytes> = vector.into();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&frozen).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&frozen));
        let mutable: FastqRecord<BytesMut> = frozen.into();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&mutable).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&mutable));
        let vector: FastqRecord<Vec<u8>> = mutable.into();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&vector).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&vector));

        let mutable: FastqRecord<BytesMut> = vector.into();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&mutable).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&mutable));
        let frozen: FastqRecord<Bytes> = mutable.into();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&frozen).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&frozen));
        let vector: FastqRecord<Vec<u8>> = frozen.into();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&vector).into_parts().0.as_ptr(),
            pointer
        );
        check_record(vector);
    }
}

#[test]
fn borrowed_conversions_copy_local_storage_and_reuse_static_storage() {
    let (vector, boxed, mutable) = {
        let data = DATA.to_vec();
        let pointer = data.as_ptr();
        let record = FastqRecord::try_from((data.as_slice(), ranges())).unwrap();
        let vector = FastqRecord::<Vec<u8>>::from(record.clone());
        let boxed = FastqRecord::<Box<[u8]>>::from(record.clone());
        let mutable = FastqRecord::<BytesMut>::from(record);
        assert_ne!(
            FastqRecord::<&[u8]>::from(&vector).into_parts().0.as_ptr(),
            pointer
        );
        assert_ne!(
            FastqRecord::<&[u8]>::from(&boxed).into_parts().0.as_ptr(),
            pointer
        );
        assert_ne!(
            FastqRecord::<&[u8]>::from(&mutable).into_parts().0.as_ptr(),
            pointer
        );
        (vector, boxed, mutable)
    };
    check_record(vector);
    check_record(boxed);
    check_record(mutable);

    let (vector, mutable, text) = {
        let data = String::from_utf8(DATA.to_vec()).unwrap();
        let pointer = data.as_ptr();
        let record = FastqRecord::try_from((data.as_str(), ranges())).unwrap();
        let vector = FastqRecord::<Vec<u8>>::from(record.clone());
        let mutable = FastqRecord::<BytesMut>::from(record.clone());
        let text = FastqRecord::<String>::from(record);
        assert_ne!(
            FastqRecord::<&[u8]>::from(&vector).into_parts().0.as_ptr(),
            pointer
        );
        assert_ne!(
            FastqRecord::<&[u8]>::from(&mutable).into_parts().0.as_ptr(),
            pointer
        );
        assert_ne!(text.as_ref().into_parts().0.as_ptr(), pointer);
        (vector, mutable, text)
    };
    check_record(vector);
    check_record(mutable);
    check_record(text);

    let bytes = FastqRecord::<Bytes>::from(FastqRecord::try_from((DATA, ranges())).unwrap());
    assert_eq!(
        FastqRecord::<&[u8]>::from(&bytes).into_parts().0.as_ptr(),
        DATA.as_ptr()
    );
    check_record(bytes);
    let text = std::str::from_utf8(DATA).unwrap();
    let bytes = FastqRecord::<Bytes>::from(FastqRecord::try_from((text, ranges())).unwrap());
    assert_eq!(
        FastqRecord::<&[u8]>::from(&bytes).into_parts().0.as_ptr(),
        text.as_ptr()
    );
    check_record(bytes);
}

#[test]
fn boxed_conversions_preserve_layout_and_reuse_unique_storage() {
    for spare_capacity in [0, 64] {
        let mut data = Vec::with_capacity(DATA.len() + spare_capacity);
        data.extend_from_slice(DATA);
        let vector = FastqRecord::try_from((data, ranges())).unwrap();
        let boxed = FastqRecord::<Box<[u8]>>::from(vector);
        let pointer = FastqRecord::<&[u8]>::from(&boxed).into_parts().0.as_ptr();
        check_record(FastqRecord::<&[u8]>::from(&boxed));
        let vector = FastqRecord::<Vec<u8>>::from(boxed);
        assert_eq!(
            FastqRecord::<&[u8]>::from(&vector).into_parts().0.as_ptr(),
            pointer
        );

        let boxed = FastqRecord::<Box<[u8]>>::from(vector);
        let bytes = FastqRecord::<Bytes>::from(boxed);
        assert_eq!(
            FastqRecord::<&[u8]>::from(&bytes).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&bytes));
        let boxed = FastqRecord::<Box<[u8]>>::from(bytes);
        let mutable = FastqRecord::<BytesMut>::from(boxed);
        assert_eq!(
            FastqRecord::<&[u8]>::from(&mutable).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&mutable));
        let boxed = FastqRecord::<Box<[u8]>>::from(mutable);
        check_record(boxed);
    }

    let text = String::from_utf8(DATA.to_vec()).unwrap();
    let record = FastqRecord::try_from((text, ranges())).unwrap();
    check_record(FastqRecord::<Box<[u8]>>::from(record));

    // Offset, shared buffers must retain bytes outside the fields as well.
    let input = [b"prefix".as_slice(), DATA].concat();
    let mut source = BytesMut::from(input.as_slice());
    let retained = source.split_to(6);
    let bytes = FastqRecord::try_from((source, ranges())).unwrap();
    check_record(FastqRecord::<Box<[u8]>>::from(bytes));
    assert_eq!(retained.as_ref(), b"prefix");

    let source = Bytes::from_static(DATA);
    let record = FastqRecord::try_from((source.clone(), ranges())).unwrap();
    let boxed = FastqRecord::<Box<[u8]>>::from(record);
    assert_ne!(
        FastqRecord::<&[u8]>::from(&boxed).into_parts().0.as_ptr(),
        source.as_ptr()
    );
    check_record(boxed);
    assert_eq!(source.as_ref(), DATA);
}

#[test]
fn string_conversions_round_trip_all_owned_byte_storage() {
    fn check<T: StringBytes + AsRef<[u8]>>()
    where
        FastqRecord<T>: From<FastqRecord<String>>,
        FastqRecord<String>: TryFrom<FastqRecord<T>, Error = FastqRecordError>,
    {
        let data = String::from_utf8(DATA.to_vec()).unwrap();
        let pointer = data.as_ptr();
        let text = <FastqRecord<String> as TryFrom<(String, FastqParts<Range<usize>>)>>::try_from(
            (data, ranges()),
        )
        .unwrap();
        let bytes = FastqRecord::<T>::from(text);
        assert_eq!(
            FastqRecord::<&[u8]>::from(&bytes).into_parts().0.as_ptr(),
            pointer
        );
        check_record(FastqRecord::<&[u8]>::from(&bytes));
        let text = FastqRecord::<String>::try_from(bytes).unwrap();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&text).into_parts().0.as_ptr(),
            pointer
        );
        check_record(text);
    }

    check::<Vec<u8>>();
    check::<Box<[u8]>>();
    check::<Bytes>();
    check::<BytesMut>();
}

#[test]
fn conversion_to_string_rejects_invalid_utf8_inside_and_outside_fields() {
    fn check<T: StringBytes + AsRef<[u8]>>(data: T, expected: &[u8], invalid_offset: usize)
    where
        FastqRecord<String>: TryFrom<FastqRecord<T>, Error = FastqRecordError>,
    {
        let record = FastqRecord::<T>::try_from((data, ranges())).unwrap();
        let FastqRecordError::InvalidUtf8(error) =
            FastqRecord::<String>::try_from(record).unwrap_err()
        else {
            panic!("expected an owned UTF-8 conversion error");
        };
        assert_eq!(error.utf8_error().valid_up_to(), invalid_offset);
        assert_eq!(error.as_bytes(), expected);
    }

    for invalid_offset in [ranges().sequence.start, DATA.len() - 1] {
        let mut invalid = DATA.to_vec();
        invalid[invalid_offset] = 255;
        check(invalid.as_slice(), &invalid, invalid_offset);
        check(invalid.clone(), &invalid, invalid_offset);
        check(invalid.clone().into_boxed_slice(), &invalid, invalid_offset);
        check(Bytes::from(invalid.clone()), &invalid, invalid_offset);
        check(BytesMut::from(invalid.as_slice()), &invalid, invalid_offset);
    }
}

#[test]
fn borrowed_text_conversions_preserve_storage_ranges_and_slice_lifetimes() {
    fn check<T: StringBytes + AsRef<[u8]>>(data: T) {
        let pointer = data.as_ref().as_ptr();
        let record = FastqRecord::try_from((data, ranges())).unwrap();
        let text: FastqRecord<&str> = (&record).try_into().unwrap();
        assert_eq!(text.clone().into_parts().0.as_ptr(), pointer);
        check_record(text);

        let text: FastqRecord<&str> = {
            let view = FastqRecord::<&[u8]>::from(&record);
            view.try_into().unwrap()
        };
        let (text, returned_ranges) = text.into_parts();
        assert_eq!(text.as_ptr(), pointer);
        assert_eq!(text.as_bytes(), DATA);
        assert_eq!(returned_ranges, ranges());
        check_record(record);
    }

    check(DATA);
    check(DATA.to_vec());
    check(Box::<[u8]>::from(DATA));
    check(Bytes::from_static(DATA));
    check(BytesMut::from(DATA));
    check(std::str::from_utf8(DATA).unwrap());
    check(String::from_utf8(DATA.to_vec()).unwrap());
    let custom = Custom(DATA.to_vec());
    check(custom.as_ref());
}

#[test]
fn borrowed_text_conversions_reject_invalid_bytes_inside_and_outside_fields() {
    for invalid_offset in [ranges().sequence.start, DATA.len() - 1] {
        let mut invalid = DATA.to_vec();
        invalid[invalid_offset] = 255;
        let record = FastqRecord::try_from((invalid, ranges())).unwrap();
        for error in [
            FastqRecord::<&str>::try_from(&record).unwrap_err(),
            FastqRecord::<&str>::try_from(FastqRecord::<&[u8]>::from(&record)).unwrap_err(),
        ] {
            let FastqRecordError::InvalidUtf8Slice(source) = error else {
                panic!("expected a borrowed UTF-8 conversion error");
            };
            assert_eq!(source.valid_up_to(), invalid_offset);
            assert_eq!(source.error_len(), Some(1));
        }
        let (data, returned_ranges) = record.into_parts();
        assert_eq!(data[invalid_offset], 255);
        assert_eq!(returned_ranges, ranges());
    }
}

#[test]
fn borrowed_storage_converts_to_a_string_that_outlives_the_source() {
    let text = {
        let data = DATA.to_vec();
        let record = FastqRecord::try_from((data.as_slice(), ranges())).unwrap();
        let text = FastqRecord::<String>::try_from(record).unwrap();
        assert_ne!(text.as_ref().into_parts().0.as_ptr(), data.as_ptr());
        text
    };
    check_record(text);
}

#[test]
fn conversion_to_string_preserves_multibyte_headers_and_byte_ranges() {
    let original = FastqRecord::<Bytes>::try_from(FastqParts {
        id: "读1".to_owned(),
        description: Some("说明".to_owned()),
        sequence: "AC".to_owned(),
        separator: "".to_owned(),
        quality: "12".to_owned(),
    })
    .unwrap();
    let expected = {
        let mut encoded = Vec::new();
        original.encode_into_vec(&mut encoded);
        encoded
    };
    let (_, expected_ranges) = FastqRecord::<&[u8]>::from(&original).into_parts();
    let text = FastqRecord::<String>::try_from(original).unwrap();
    assert_eq!(text.id(), "读1".as_bytes());
    assert_eq!(text.description(), Some("说明".as_bytes()));
    assert_eq!(text.sequence(), b"AC");
    assert_eq!(text.quality(), b"12");
    assert_eq!(
        {
            let mut encoded = Vec::new();
            text.encode_into_vec(&mut encoded);
            encoded
        },
        expected
    );
    assert_eq!(text.into_parts().1, expected_ranges);
}

#[test]
fn from_bytes_copies_shared_and_static_storage() {
    for source in [Bytes::copy_from_slice(DATA), Bytes::from_static(DATA)] {
        let frozen = FastqRecord::<Bytes>::try_from((source.clone(), ranges())).unwrap();
        let vector = FastqRecord::<Vec<u8>>::from(frozen.clone());
        assert_ne!(
            FastqRecord::<&[u8]>::from(&vector).into_parts().0.as_ptr(),
            source.as_ptr()
        );
        check_record(FastqRecord::<&[u8]>::from(&vector));
        let (mut data, ranges) = vector.into_parts();
        data[ranges.sequence.start] = b'T';
        assert_eq!(source.as_ref(), DATA);

        let mutable = FastqRecord::<BytesMut>::from(frozen);
        assert_ne!(
            FastqRecord::<&[u8]>::from(&mutable).into_parts().0.as_ptr(),
            source.as_ptr()
        );
        check_record(FastqRecord::<&[u8]>::from(&mutable));

        let (mut data, ranges) = mutable.into_parts();
        data[ranges.sequence.start] = b'T';
        assert_eq!(source.as_ref(), DATA);
    }
}

#[test]
fn converting_offset_storage_to_vectors_preserves_fields_and_shared_data() {
    let input = [b"prefix".as_slice(), DATA].concat();
    for shared in [false, true] {
        let mut source = Bytes::copy_from_slice(&input);
        let prefix = source.split_to(6);
        let prefix = shared.then_some(prefix);
        let record = FastqRecord::try_from((source, ranges())).unwrap();
        check_record(FastqRecord::<Vec<u8>>::from(record));
        if let Some(prefix) = prefix {
            assert_eq!(prefix.as_ref(), b"prefix");
        }

        let mut source = BytesMut::from(input.as_slice());
        let prefix = source.split_to(6);
        let prefix = shared.then_some(prefix);
        let record = FastqRecord::try_from((source, ranges())).unwrap();
        check_record(FastqRecord::<Vec<u8>>::from(record));
        if let Some(prefix) = prefix {
            assert_eq!(prefix.as_ref(), b"prefix");
        }
    }
}

#[test]
fn overlapping_ranges_and_empty_fields_are_valid() {
    let overlapping = FastqRecord::try_from((
        b"AC".as_slice(),
        FastqParts {
            id: 0..2,
            description: None,
            sequence: 0..2,
            separator: 2..2,
            quality: 0..2,
        },
    ))
    .unwrap();
    assert_eq!(
        {
            let mut encoded = Vec::new();
            overlapping.encode_into_vec(&mut encoded);
            encoded
        },
        b"@AC\nAC\n+\nAC\n"
    );
    for (description, expected) in [
        (None, b"@\n\n+\n\n".as_slice()),
        (Some(0..0), b"@ \n\n+\n\n".as_slice()),
    ] {
        let record = FastqRecord::try_from((
            b"".as_slice(),
            FastqParts {
                id: 0..0,
                description,
                sequence: 0..0,
                separator: 0..0,
                quality: 0..0,
            },
        ))
        .unwrap();
        assert_eq!(
            {
                let mut encoded = Vec::new();
                record.encode_into_vec(&mut encoded);
                encoded
            },
            expected
        );
    }
}
