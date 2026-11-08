use std::cell::Cell;
use std::convert::Infallible;
use std::error::Error;
use std::ops::Deref;
use std::str::Utf8Error;
use std::string::FromUtf8Error;

use bytes::{Bytes, BytesMut};

use biorill_fastx::{FastqField, FastqFieldError, FastqParts};
use stringbytes::StringBytes;

struct CustomBytes(Vec<u8>);

impl AsRef<[u8]> for CustomBytes {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

fn borrow_parts<T: AsRef<[u8]>>(parts: &FastqParts<T>) -> FastqParts<&[u8]> {
    FastqParts::new(
        parts.id.as_ref(),
        parts.description.as_ref().map(AsRef::as_ref),
        parts.sequence.as_ref(),
        parts.separator.as_ref(),
        parts.quality.as_ref(),
    )
}

fn parts() -> FastqParts<&'static [u8]> {
    FastqParts {
        id: b"read1",
        description: Some(b"note"),
        sequence: b"A\xff",
        separator: b"read1",
        quality: b"12",
    }
}

fn vectors(spare: usize) -> FastqField<Vec<u8>> {
    let allocate = |bytes: &[u8]| {
        let mut data = Vec::with_capacity(bytes.len() + spare);
        data.extend_from_slice(bytes);
        data
    };
    let input = parts();
    FastqField::try_from(FastqParts {
        id: allocate(input.id),
        description: input.description.map(allocate),
        sequence: allocate(input.sequence),
        separator: allocate(input.separator),
        quality: allocate(input.quality),
    })
    .unwrap()
}

fn check<T: StringBytes + AsRef<[u8]>>(fields: &FastqField<T>) -> FastqParts<*const u8> {
    let borrowed = FastqField::<&[u8]>::from(fields).into_parts();
    assert_eq!(borrowed, parts());
    FastqParts {
        id: borrowed.id.as_ptr(),
        description: borrowed.description.map(|bytes| bytes.as_ptr()),
        sequence: borrowed.sequence.as_ptr(),
        separator: borrowed.separator.as_ptr(),
        quality: borrowed.quality.as_ptr(),
    }
}

#[test]
fn as_ref_borrows_original_components_and_as_deref_follows_one_level() {
    // Borrowing the components needs no conversion or Clone bound.
    fn borrow<T>(fields: &FastqField<T>) -> FastqField<&T> {
        fields.as_ref()
    }

    fn check<T>(wrap: impl Fn(&'static str) -> T)
    where
        T: StringBytes + AsRef<[u8]> + Deref,
        T::Target: StringBytes + AsRef<[u8]>,
        FastqField<T>: TryFrom<FastqParts<T>, Error = FastqFieldError>,
    {
        for description in [None, Some(""), Some("说明")] {
            let fields = FastqField::try_from(FastqParts::new(
                wrap("read1"),
                description.map(&wrap),
                wrap("AC"),
                wrap(""),
                wrap("12"),
            ))
            .unwrap();
            let borrowed: FastqField<&T> = borrow(&fields);
            // Dereferencing &T yields T, without following T's own Deref again.
            let once: FastqField<&T> = borrowed.as_deref();
            for view in [&borrowed, &once] {
                for (original, reference) in [
                    (fields.id(), *view.id()),
                    (fields.sequence(), *view.sequence()),
                    (fields.separator(), *view.separator()),
                    (fields.quality(), *view.quality()),
                ] {
                    assert!(std::ptr::eq(original, reference));
                }
                assert_eq!(view.description().is_some(), description.is_some());
                if let Some(original) = fields.description() {
                    assert!(std::ptr::eq(original, *view.description().unwrap()));
                }
            }

            let dereferenced: FastqField<&T::Target> = fields.as_deref();
            for (original, target) in [
                (fields.id(), *dereferenced.id()),
                (fields.sequence(), *dereferenced.sequence()),
                (fields.separator(), *dereferenced.separator()),
                (fields.quality(), *dereferenced.quality()),
            ] {
                assert_eq!(original.as_ref(), target.as_ref());
                assert_eq!(original.as_ref().as_ptr(), target.as_ref().as_ptr());
            }
            assert_eq!(dereferenced.description().is_some(), description.is_some());
            if let Some(original) = fields.description() {
                let target = *dereferenced.description().unwrap();
                assert_eq!(original.as_ref(), target.as_ref());
                assert_eq!(original.as_ref().as_ptr(), target.as_ref().as_ptr());
            }
        }
    }

    check(str::as_bytes);
    check(|text| Vec::from(text.as_bytes()));
    check(|text| Box::<[u8]>::from(text.as_bytes()));
    check(Bytes::from);
    check(BytesMut::from);
    check(|text| text);
    check(String::from);
}

#[test]
fn all_owned_byte_conversions_preserve_fields_and_reuse_unique_buffers() {
    for spare in [0, 32] {
        let fields = vectors(spare);
        let pointers = check(&fields);
        let fields: FastqField<Bytes> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<BytesMut> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Vec<u8>> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Box<[u8]>> = fields.into();
        let pointers = check(&fields);
        let fields: FastqField<Vec<u8>> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Box<[u8]>> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Bytes> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Box<[u8]>> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<BytesMut> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Box<[u8]>> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Vec<u8>> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<BytesMut> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Bytes> = fields.into();
        assert_eq!(check(&fields), pointers);
        let fields: FastqField<Vec<u8>> = fields.into();
        assert_eq!(check(&fields), pointers);
    }
}

#[test]
fn shared_and_static_fields_are_copied_into_independent_mutable_storage() {
    let input = parts();
    let static_fields = FastqField::try_from(FastqParts {
        id: Bytes::from_static(input.id),
        description: input.description.map(Bytes::from_static),
        sequence: Bytes::from_static(input.sequence),
        separator: Bytes::from_static(input.separator),
        quality: Bytes::from_static(input.quality),
    })
    .unwrap();
    for fields in [FastqField::<Bytes>::from(vectors(0)), static_fields] {
        let source = check(&fields);
        let vector = FastqField::<Vec<u8>>::from(fields.clone());
        assert_ne!(check(&vector).sequence, source.sequence);
        let mut owned = vector.into_parts();
        owned.sequence[0] = b'T';
        assert_eq!(fields.sequence().as_ref(), input.sequence);

        let mutable = FastqField::<BytesMut>::from(fields.clone());
        assert_ne!(check(&mutable).sequence, source.sequence);
        let mut owned = mutable.into_parts();
        owned.sequence[0] = b'T';
        assert_eq!(fields.sequence().as_ref(), input.sequence);

        let boxed = FastqField::<Box<[u8]>>::from(fields.clone());
        assert_ne!(check(&boxed).sequence, source.sequence);
        let mut owned = boxed.into_parts();
        owned.sequence[0] = b'T';
        assert_eq!(fields.sequence().as_ref(), input.sequence);
    }

    let mut storage = BytesMut::from(b"read1noteA\xffread112tail".as_slice());
    let fields = FastqField::try_from(FastqParts {
        id: storage.split_to(5),
        description: Some(storage.split_to(4)),
        sequence: storage.split_to(2),
        separator: storage.split_to(5),
        quality: storage.split_to(2),
    })
    .unwrap();
    let source = check(&fields);
    let boxed = FastqField::<Box<[u8]>>::from(fields);
    let copied = check(&boxed);
    for (source, copied) in [
        (source.id, copied.id),
        (source.description.unwrap(), copied.description.unwrap()),
        (source.sequence, copied.sequence),
        (source.separator, copied.separator),
        (source.quality, copied.quality),
    ] {
        assert_ne!(source, copied);
    }
    // Keeping the tail alive keeps the original allocation shared during conversion.
    assert_eq!(storage.as_ref(), b"tail");
}

#[test]
fn try_into_string_round_trips_all_owned_byte_types_and_optional_descriptions() {
    fn check<T: StringBytes + AsRef<[u8]>>()
    where
        FastqField<T>: From<FastqField<String>>,
        FastqField<String>: TryFrom<FastqField<T>, Error = FastqFieldError>,
    {
        for description in [None, Some(""), Some("说明")] {
            for (sequence, quality) in [("AC", "12"), ("", "")] {
                let input = FastqParts {
                    id: String::from("读1"),
                    description: description.map(String::from),
                    sequence: String::from(sequence),
                    separator: String::from("read1"),
                    quality: String::from(quality),
                };
                let fields =
                    <FastqField<String> as TryFrom<FastqParts<String>>>::try_from(input.clone())
                        .unwrap();
                let pointer = fields.id().as_ptr();
                let bytes = FastqField::<T>::from(fields);
                assert_eq!(
                    FastqField::<&[u8]>::from(&bytes).into_parts(),
                    borrow_parts(&input)
                );
                assert_eq!(bytes.id().as_ref().as_ptr(), pointer);
                let text: FastqField<String> = bytes.try_into().unwrap();
                assert_eq!(text.id().as_ptr(), pointer);
                assert_eq!(text.into_parts(), input);
            }
        }
    }

    check::<Vec<u8>>();
    check::<Bytes>();
    check::<BytesMut>();
    check::<Box<[u8]>>();
}

#[test]
fn boxed_fields_convert_to_string_without_copying_their_bytes() {
    for description in [None, Some(b"".as_slice()), Some("说明".as_bytes())] {
        for (sequence, quality) in [(b"AC".as_slice(), b"12".as_slice()), (b"", b"")] {
            let expected = FastqParts {
                id: "读1".as_bytes(),
                description,
                sequence,
                separator: b"read1".as_slice(),
                quality,
            };
            let fields = FastqField::try_from(FastqParts {
                id: Box::<[u8]>::from(expected.id),
                description: expected.description.map(Box::<[u8]>::from),
                sequence: Box::<[u8]>::from(expected.sequence),
                separator: Box::<[u8]>::from(expected.separator),
                quality: Box::<[u8]>::from(expected.quality),
            })
            .unwrap();
            let pointers = [
                fields.id().as_ptr(),
                fields.sequence().as_ptr(),
                fields.separator().as_ptr(),
                fields.quality().as_ptr(),
            ];
            let description_pointer = fields.description().map(|bytes| bytes.as_ptr());
            let text: FastqField<String> = fields.try_into().unwrap();
            assert_eq!(FastqField::<&[u8]>::from(&text).into_parts(), expected);
            assert_eq!(
                [
                    text.id().as_ptr(),
                    text.sequence().as_ptr(),
                    text.separator().as_ptr(),
                    text.quality().as_ptr(),
                ],
                pointers
            );
            assert_eq!(
                text.description().map(|text| text.as_ptr()),
                description_pointer
            );
        }
    }
}

#[test]
fn string_fields_borrow_as_str_without_copying() {
    for description in [None, Some(""), Some("说明")] {
        for (sequence, quality) in [("AC", "12"), ("", "")] {
            let expected = FastqParts {
                id: String::from("读1"),
                description: description.map(String::from),
                sequence: String::from(sequence),
                separator: String::from("read1"),
                quality: String::from(quality),
            };
            let source = FastqField::try_from(expected.clone()).unwrap();
            let borrowed = FastqField::<&str>::from(&source);
            let via_into: FastqField<&str> = (&source).into();
            let via_try_into: Result<FastqField<&str>, Infallible> = (&source).try_into();
            assert_eq!(borrowed, via_into);
            assert_eq!(borrowed, via_try_into.unwrap());
            assert_eq!(borrowed, source.as_deref());
            let reborrowed = borrowed.as_deref();
            assert_eq!(reborrowed, borrowed);
            for fields in [&borrowed, &reborrowed] {
                assert_eq!(fields.id().as_ptr(), source.id().as_ptr());
                assert_eq!(
                    fields.description().map(|text| text.as_ptr()),
                    source.description().map(|text| text.as_ptr())
                );
                assert_eq!(fields.sequence().as_ptr(), source.sequence().as_ptr());
                assert_eq!(fields.separator().as_ptr(), source.separator().as_ptr());
                assert_eq!(fields.quality().as_ptr(), source.quality().as_ptr());
            }
            assert_eq!(
                FastqField::<&[u8]>::from(&borrowed).into_parts(),
                borrow_parts(&expected)
            );
        }
    }
}

#[test]
fn try_into_str_borrows_byte_storage_and_preserves_field_contents() {
    fn check<T: StringBytes + AsRef<[u8]>>(wrap: impl Fn(&str) -> T)
    where
        FastqField<T>: TryFrom<FastqParts<T>, Error = FastqFieldError>,
        for<'a> FastqField<&'a str>: TryFrom<&'a FastqField<T>, Error = FastqFieldError>,
    {
        for description in [None, Some(""), Some("说明")] {
            for (sequence, quality) in [("AC", "12"), ("", "")] {
                let bytes = FastqField::<T>::try_from(FastqParts {
                    id: wrap("读1"),
                    description: description.map(&wrap),
                    sequence: wrap(sequence),
                    separator: wrap("read1"),
                    quality: wrap(quality),
                })
                .unwrap();
                let borrowed: FastqField<&str> = (&bytes).try_into().unwrap();
                let consumed: FastqField<&str> = {
                    let slices = FastqField::<&[u8]>::from(&bytes);
                    slices.try_into().unwrap()
                };
                assert_eq!(borrowed, consumed);
                let slices = FastqField::<&[u8]>::from(&bytes);
                let from_slices: FastqField<&str> = (&slices).try_into().unwrap();
                assert_eq!(borrowed, from_slices);
                assert_eq!(
                    FastqField::<&[u8]>::from(&borrowed).into_parts(),
                    slices.into_parts()
                );
                for text in [&borrowed, &consumed] {
                    for (text, data) in [
                        (*text.id(), bytes.id().as_ref()),
                        (*text.sequence(), bytes.sequence().as_ref()),
                        (*text.separator(), bytes.separator().as_ref()),
                        (*text.quality(), bytes.quality().as_ref()),
                    ] {
                        assert_eq!(text.as_ptr(), data.as_ptr());
                    }
                    assert_eq!(
                        text.description().map(|text| text.as_ptr()),
                        bytes.description().map(|data| data.as_ref().as_ptr())
                    );
                }
            }
        }
    }

    check(|text| text.as_bytes().to_vec());
    check(|text| Bytes::copy_from_slice(text.as_bytes()));
    check(|text| BytesMut::from(text));
    check(|text| text.as_bytes().to_vec().into_boxed_slice());
}

#[test]
fn try_into_str_reports_each_invalid_field_without_changing_the_source() {
    fn check<T: StringBytes + AsRef<[u8]>>(fields: FastqField<T>, name: &str, valid_prefix: usize)
    where
        for<'a> FastqField<&'a str>: TryFrom<&'a FastqField<T>, Error = FastqFieldError>,
    {
        let slices = FastqField::<&[u8]>::from(&fields);
        let before = slices.clone().into_parts();
        let results: [Result<FastqField<&str>, FastqFieldError>; 3] = [
            (&slices).try_into(),
            (&fields).try_into(),
            FastqField::<&[u8]>::from(&fields).try_into(),
        ];
        for result in results {
            let error = result.unwrap_err();
            let FastqFieldError::InvalidUtf8Slice { field, source } = &error else {
                panic!("expected a borrowed UTF-8 error");
            };
            assert_eq!(*field, name);
            assert_eq!(source.valid_up_to(), valid_prefix);
            assert_eq!(source.error_len(), Some(1));
            assert_eq!(
                error.source().unwrap().downcast_ref::<Utf8Error>(),
                Some(source)
            );
        }
        assert_eq!(FastqField::<&[u8]>::from(&fields).into_parts(), before);
    }

    for (name, invalid, valid_prefix) in [
        ("identifier", b"\xff".as_slice(), 0),
        ("description", b"\xff".as_slice(), 0),
        ("sequence", b"A\xff".as_slice(), 1),
        ("separator", b"\xff".as_slice(), 0),
        ("quality", b"1\xff".as_slice(), 1),
    ] {
        let mut input = FastqParts {
            sequence: b"AC".as_slice(),
            ..parts()
        };
        match name {
            "identifier" => input.id = invalid,
            "description" => input.description = Some(invalid),
            "sequence" => input.sequence = invalid,
            "separator" => input.separator = invalid,
            "quality" => input.quality = invalid,
            _ => unreachable!(),
        }
        let fields = FastqField::try_from(FastqParts {
            id: input.id.to_vec(),
            description: input.description.map(Vec::from),
            sequence: input.sequence.to_vec(),
            separator: input.separator.to_vec(),
            quality: input.quality.to_vec(),
        })
        .unwrap();
        check(fields.clone(), name, valid_prefix);
        check(
            FastqField::<Bytes>::from(fields.clone()),
            name,
            valid_prefix,
        );
        check(FastqField::<BytesMut>::from(fields), name, valid_prefix);
        let custom = FastqParts::new(
            CustomBytes(input.id.to_vec()),
            input.description.map(|bytes| CustomBytes(bytes.to_vec())),
            CustomBytes(input.sequence.to_vec()),
            CustomBytes(input.separator.to_vec()),
            CustomBytes(input.quality.to_vec()),
        );
        let borrowed = FastqField::try_from(borrow_parts(&custom)).unwrap();
        check(borrowed, name, valid_prefix);
    }
}

#[test]
fn try_into_string_copies_nonstatic_slices_and_outlives_the_source() {
    for description in [None, Some(""), Some("说明")] {
        for (sequence, quality) in [("AC", "12"), ("", ""), ("é", "12")] {
            let expected = FastqParts {
                id: String::from("读1"),
                description: description.map(String::from),
                sequence: String::from(sequence),
                separator: String::from("read1"),
                quality: String::from(quality),
            };
            let owned = {
                let source = FastqField::try_from(expected.clone()).unwrap();
                let slices = FastqField::<&[u8]>::from(&source);
                let owned: FastqField<String> = slices.try_into().unwrap();
                assert_eq!(
                    FastqField::<&[u8]>::from(&owned),
                    FastqField::<&[u8]>::from(&source)
                );
                for (copied, original) in [
                    (owned.id(), source.id()),
                    (owned.sequence(), source.sequence()),
                    (owned.separator(), source.separator()),
                    (owned.quality(), source.quality()),
                ] {
                    if !original.is_empty() {
                        assert_ne!(copied.as_ptr(), original.as_ptr());
                    }
                }
                if let Some(original) = source.description().filter(|text| !text.is_empty()) {
                    assert_ne!(owned.description().unwrap().as_ptr(), original.as_ptr());
                }
                owned
            };
            assert_eq!(owned.into_parts(), expected);
        }
    }
}

#[test]
fn native_str_conversions_support_infallible_try_into_with_nonstatic_fields() {
    fn check<T: StringBytes + AsRef<[u8]>>()
    where
        for<'a> FastqField<T>: From<FastqField<&'a str>>,
    {
        for description in [None, Some(""), Some("说明")] {
            for (sequence, quality) in [("AC", "12"), ("", "")] {
                let expected = FastqParts {
                    id: String::from("读1"),
                    description: description.map(String::from),
                    sequence: String::from(sequence),
                    separator: String::from("read1"),
                    quality: String::from(quality),
                };
                let owned = {
                    let source = FastqField::try_from(expected.clone()).unwrap();
                    let result: Result<FastqField<T>, Infallible> = source.as_deref().try_into();
                    result.unwrap()
                };
                assert_eq!(
                    FastqField::<&[u8]>::from(&owned).into_parts(),
                    borrow_parts(&expected)
                );
            }
        }
    }

    check::<String>();
    check::<Vec<u8>>();
    check::<BytesMut>();
}

#[test]
fn native_from_static_str_into_bytes_preserves_field_pointers() {
    for description in [None, Some(""), Some("说明")] {
        for (sequence, quality) in [("AC", "12"), ("", "")] {
            let input = FastqParts {
                id: "读1",
                description,
                sequence,
                separator: "read1",
                quality,
            };
            let source = FastqField::<&'static str>::try_from(input.clone()).unwrap();
            let fields = FastqField::<Bytes>::from(source);
            let borrowed = FastqField::<&[u8]>::from(&fields).into_parts();
            assert_eq!(borrowed, borrow_parts(&input));
            for (source, target) in [
                (input.id, borrowed.id),
                (input.sequence, borrowed.sequence),
                (input.separator, borrowed.separator),
                (input.quality, borrowed.quality),
            ] {
                if !source.is_empty() {
                    assert_eq!(source.as_ptr(), target.as_ptr());
                }
            }
            if let Some(description) = input.description.filter(|text| !text.is_empty()) {
                assert_eq!(description.as_ptr(), borrowed.description.unwrap().as_ptr());
            }
        }
    }
}

#[test]
fn native_slice_conversions_copy_nonstatic_fields_without_requiring_utf8() {
    fn check<T: StringBytes + AsRef<[u8]>>()
    where
        for<'a> FastqField<T>: From<FastqField<&'a [u8]>>,
    {
        for description in [None, Some(b"".as_slice()), Some(b"\xff".as_slice())] {
            for (sequence, quality) in [(b"A\xff".as_slice(), b"12".as_slice()), (b"", b"")] {
                let expected = FastqParts {
                    description,
                    sequence,
                    quality,
                    ..parts()
                };
                let owned = {
                    let source = FastqField::try_from(FastqParts {
                        id: expected.id.to_vec(),
                        description: expected.description.map(Vec::from),
                        sequence: expected.sequence.to_vec(),
                        separator: expected.separator.to_vec(),
                        quality: expected.quality.to_vec(),
                    })
                    .unwrap();
                    let result: Result<FastqField<T>, Infallible> =
                        FastqField::<&[u8]>::from(&source).try_into();
                    let owned = result.unwrap();
                    assert_ne!(owned.id().as_ref().as_ptr(), source.id().as_ptr());
                    owned
                };
                assert_eq!(FastqField::<&[u8]>::from(&owned).into_parts(), expected);
            }
        }
    }

    check::<Vec<u8>>();
    check::<BytesMut>();
    check::<Box<[u8]>>();
}

#[test]
fn native_from_static_slices_into_bytes_preserves_field_pointers() {
    for description in [None, Some(b"".as_slice()), Some(b"\xff".as_slice())] {
        for (sequence, quality) in [(b"A\xff".as_slice(), b"12".as_slice()), (b"", b"")] {
            let input = FastqParts {
                description,
                sequence,
                quality,
                ..parts()
            };
            let source = FastqField::<&'static [u8]>::try_from(input.clone()).unwrap();
            let fields = FastqField::<Bytes>::from(source);
            let borrowed = FastqField::<&[u8]>::from(&fields).into_parts();
            assert_eq!(borrowed, input);
            for (source, target) in [
                (input.id, borrowed.id),
                (input.sequence, borrowed.sequence),
                (input.separator, borrowed.separator),
                (input.quality, borrowed.quality),
            ] {
                if !source.is_empty() {
                    assert_eq!(source.as_ptr(), target.as_ptr());
                }
            }
            if let Some(description) = input.description.filter(|bytes| !bytes.is_empty()) {
                assert_eq!(description.as_ptr(), borrowed.description.unwrap().as_ptr());
            }
        }
    }
}

#[test]
fn borrowed_fields_do_not_revisit_custom_storage() {
    struct View {
        text: String,
        byte_reads: Cell<usize>,
    }

    impl AsRef<[u8]> for View {
        fn as_ref(&self) -> &[u8] {
            self.byte_reads.set(self.byte_reads.get() + 1);
            self.text.as_bytes()
        }
    }

    let view = |text: &str| View {
        text: text.to_owned(),
        byte_reads: Cell::new(0),
    };
    let parts = FastqParts::new(
        view("read1"),
        Some(view("note")),
        view("AC"),
        view("read1"),
        view("12"),
    );
    let fields = FastqField::try_from(borrow_parts(&parts)).unwrap();
    let borrowed = FastqField::<&[u8]>::from(&fields);
    let text: FastqField<&str> = (&fields).try_into().unwrap();
    assert_eq!(*text.id(), "read1");
    assert_eq!(text.id().as_ptr(), parts.id.text.as_ptr());
    for (source, borrowed) in [
        (&parts.id, *borrowed.id()),
        (
            parts.description.as_ref().unwrap(),
            *borrowed.description().unwrap(),
        ),
        (&parts.sequence, *borrowed.sequence()),
        (&parts.separator, *borrowed.separator()),
        (&parts.quality, *borrowed.quality()),
    ] {
        assert_eq!(source.text.as_bytes(), borrowed);
        assert_eq!(source.text.as_ptr(), borrowed.as_ptr());
        assert_eq!(source.byte_reads.get(), 1);
    }
    assert_eq!(
        FastqField::<&[u8]>::from(&fields).into_parts(),
        FastqField::<&[u8]>::from(&text).into_parts()
    );
    assert_eq!(parts.id.byte_reads.get(), 1);
}

#[test]
fn try_into_string_identifies_each_invalid_field_and_preserves_its_bytes() {
    fn check<T>(fields: FastqField<T>, name: &str, bytes: &[u8], valid_prefix: usize)
    where
        FastqField<String>: TryFrom<FastqField<T>, Error = FastqFieldError>,
    {
        let result: Result<FastqField<String>, FastqFieldError> = fields.try_into();
        let error = result.unwrap_err();
        let FastqFieldError::InvalidUtf8 { field, source } = &error else {
            panic!("expected a UTF-8 conversion error");
        };
        assert_eq!(*field, name);
        assert_eq!(source.as_bytes(), bytes);
        assert_eq!(source.utf8_error().valid_up_to(), valid_prefix);
        assert_eq!(source.utf8_error().error_len(), Some(1));
        assert_eq!(
            error.source().unwrap().downcast_ref::<FromUtf8Error>(),
            Some(source)
        );
    }

    for (name, invalid, valid_prefix) in [
        ("identifier", b"\xff".as_slice(), 0),
        ("description", b"\xff".as_slice(), 0),
        ("sequence", b"A\xff".as_slice(), 1),
        ("separator", b"\xff".as_slice(), 0),
        ("quality", b"1\xff".as_slice(), 1),
    ] {
        let mut input = FastqParts {
            sequence: b"AC".as_slice(),
            ..parts()
        };
        match name {
            "identifier" => input.id = invalid,
            "description" => input.description = Some(invalid),
            "sequence" => input.sequence = invalid,
            "separator" => input.separator = invalid,
            "quality" => input.quality = invalid,
            _ => unreachable!(),
        }
        let fields = FastqField::try_from(FastqParts {
            id: input.id.to_vec(),
            description: input.description.map(Vec::from),
            sequence: input.sequence.to_vec(),
            separator: input.separator.to_vec(),
            quality: input.quality.to_vec(),
        })
        .unwrap();
        check(
            FastqField::<&[u8]>::from(&fields),
            name,
            invalid,
            valid_prefix,
        );
        let boxed = FastqField::try_from(FastqParts {
            id: Box::<[u8]>::from(input.id),
            description: input.description.map(Box::<[u8]>::from),
            sequence: Box::<[u8]>::from(input.sequence),
            separator: Box::<[u8]>::from(input.separator),
            quality: Box::<[u8]>::from(input.quality),
        })
        .unwrap();
        check(boxed, name, invalid, valid_prefix);
        check(fields.clone(), name, invalid, valid_prefix);
        check(
            FastqField::<Bytes>::from(fields.clone()),
            name,
            invalid,
            valid_prefix,
        );
        check(
            FastqField::<BytesMut>::from(fields),
            name,
            invalid,
            valid_prefix,
        );
    }
}
