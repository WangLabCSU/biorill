use std::borrow::Cow;
use std::cell::Cell;
use std::ops::Deref;
use std::rc::Rc;
use std::sync::Arc;

use bytes::{Bytes, BytesMut};

use biorill_fastx::{FastqField, FastqFieldError, FastqParts, FastqRecord, LengthMismatch};
use stringbytes::StringBytes;

fn parts() -> FastqParts<&'static [u8]> {
    FastqParts {
        id: b"read1",
        description: Some(b"description"),
        sequence: b"AC",
        separator: b"read1",
        quality: b"12",
    }
}

#[test]
fn parts_as_ref_borrows_components_without_calling_their_traits() {
    fn borrow<T>(parts: &FastqParts<T>) -> FastqParts<&T> {
        parts.as_ref()
    }

    struct Component {
        value: i32,
        derefs: Cell<usize>,
    }

    impl Deref for Component {
        type Target = i32;

        fn deref(&self) -> &i32 {
            self.derefs.set(self.derefs.get() + 1);
            &self.value
        }
    }

    let component = |value| Component {
        value,
        derefs: Cell::new(0),
    };
    for description in [None, Some(5)] {
        let parts = FastqParts::new(
            component(1),
            description.map(component),
            component(2),
            component(3),
            component(4),
        );
        let borrowed: FastqParts<&Component> = borrow(&parts);
        let once: FastqParts<&Component> = borrowed.as_deref();
        for view in [&borrowed, &once] {
            for (original, reference) in [
                (&parts.id, view.id),
                (&parts.sequence, view.sequence),
                (&parts.separator, view.separator),
                (&parts.quality, view.quality),
            ] {
                assert!(std::ptr::eq(original, reference));
                assert_eq!(original.derefs.get(), 0);
            }
            assert_eq!(view.description.is_some(), description.is_some());
            if let Some(original) = &parts.description {
                assert!(std::ptr::eq(original, view.description.unwrap()));
                assert_eq!(original.derefs.get(), 0);
            }
        }

        let values: FastqParts<&i32> = parts.as_deref();
        for (original, value) in [
            (&parts.id, values.id),
            (&parts.sequence, values.sequence),
            (&parts.separator, values.separator),
            (&parts.quality, values.quality),
        ] {
            assert!(std::ptr::eq(&original.value, value));
            assert_eq!(original.derefs.get(), 1);
        }
        assert_eq!(values.description.copied(), description);
        if let Some(original) = &parts.description {
            assert!(std::ptr::eq(&original.value, values.description.unwrap()));
            assert_eq!(original.derefs.get(), 1);
        }
    }
}

#[test]
fn parts_as_deref_borrows_text_bytes_and_sized_targets_without_validation() {
    fn check<T: Deref>(wrap: impl Fn(&'static str) -> T)
    where
        T::Target: AsRef<[u8]>,
    {
        for description in [None, Some(""), Some("note\nextra")] {
            let parts = FastqParts::new(
                wrap("read 1"),
                description.map(&wrap),
                wrap("AC"),
                wrap(""),
                wrap("1"),
            );
            let view = parts.as_deref();
            assert_eq!(view.id.as_ref(), b"read 1");
            assert_eq!(view.sequence.as_ref(), b"AC");
            assert_eq!(view.quality.as_ref(), b"1");
            for (original, target) in [
                (&parts.id, view.id),
                (&parts.sequence, view.sequence),
                (&parts.separator, view.separator),
                (&parts.quality, view.quality),
            ] {
                assert!(std::ptr::eq(original.deref(), target));
            }
            assert_eq!(
                view.description.map(AsRef::as_ref),
                description.map(str::as_bytes)
            );
            if let Some(original) = &parts.description {
                assert!(std::ptr::eq(original.deref(), view.description.unwrap()));
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

    let numbers = FastqParts::new(Box::new(1), None, Box::new(2), Box::new(3), Box::new(4));
    let borrowed: FastqParts<&i32> = numbers.as_deref();
    assert!(std::ptr::eq(borrowed.id, numbers.id.as_ref()));
    assert_eq!(borrowed.description, None);
    assert_eq!(*borrowed.sequence, 2);
    assert_eq!(*borrowed.separator, 3);
    assert_eq!(*borrowed.quality, 4);
}

#[test]
fn parts_conversion_accepts_common_storage_and_borrows_without_copying() {
    fn check<T>(wrap: impl Fn(&'static str) -> T)
    where
        T: Deref,
        T::Target: StringBytes + AsRef<[u8]>,
        FastqField<T>: TryFrom<FastqParts<T>, Error = FastqFieldError>,
        for<'a> FastqField<&'a T>: TryFrom<FastqParts<&'a T>, Error = FastqFieldError>,
    {
        for description in [None, Some(""), Some("说明")] {
            let input = FastqParts::new(
                wrap("read1"),
                description.map(&wrap),
                wrap("AC"),
                wrap(""),
                wrap("12"),
            );
            let source = input.as_deref();
            let pointers = [
                source.id.as_ref().as_ptr(),
                source.sequence.as_ref().as_ptr(),
                source.separator.as_ref().as_ptr(),
                source.quality.as_ref().as_ptr(),
            ];
            let description_pointer = source.description.map(|value| value.as_ref().as_ptr());
            let borrowed: FastqField<&T> = FastqField::try_from(input.as_ref()).unwrap();
            for (original, reference) in [
                (&input.id, *borrowed.id()),
                (&input.sequence, *borrowed.sequence()),
                (&input.separator, *borrowed.separator()),
                (&input.quality, *borrowed.quality()),
            ] {
                assert!(std::ptr::eq(original, reference));
            }
            if let Some(original) = &input.description {
                assert!(std::ptr::eq(original, *borrowed.description().unwrap()));
            } else {
                assert!(borrowed.description().is_none());
            }

            let fields = FastqField::try_from(input).unwrap();
            let view = fields.as_deref();
            let components = view.into_parts();
            assert_eq!(
                pointers,
                [
                    components.id.as_ref().as_ptr(),
                    components.sequence.as_ref().as_ptr(),
                    components.separator.as_ref().as_ptr(),
                    components.quality.as_ref().as_ptr(),
                ]
            );
            assert_eq!(
                description_pointer,
                components.description.map(|value| value.as_ref().as_ptr())
            );
            let mut encoded = Vec::new();
            fields.as_deref().encode_into_vec(&mut encoded);
            let header =
                description.map_or(String::from("@read1"), |text| format!("@read1 {text}"));
            assert_eq!(encoded, format!("{header}\nAC\n+\n12\n").as_bytes());
        }

        for (id, sequence, quality, expected) in [
            ("read1", "A\n", "12", "line feed"),
            ("read 1", "AC", "12", "identifier"),
            ("read1", "AC", "1", "length"),
        ] {
            let input = FastqParts::new(wrap(id), None, wrap(sequence), wrap(""), wrap(quality));
            let borrowed_error = FastqField::try_from(input.as_ref()).err().unwrap();
            let owned_error = FastqField::try_from(input).err().unwrap();
            for error in [borrowed_error, owned_error] {
                match expected {
                    "line feed" => assert!(matches!(
                        error,
                        FastqFieldError::LineFeed { field: "sequence" }
                    )),
                    "identifier" => assert!(matches!(error, FastqFieldError::InvalidId)),
                    "length" => assert!(matches!(error, FastqFieldError::LengthMismatch(_))),
                    _ => unreachable!(),
                }
            }
        }
    }

    check(str::as_bytes);
    check(|text| Vec::from(text.as_bytes()));
    check(|text| Box::<[u8]>::from(text.as_bytes()));
    check(Bytes::from);
    check(BytesMut::from);
    check(|text| Cow::Borrowed(text.as_bytes()));
    check(|text| Cow::<'static, [u8]>::Owned(Vec::from(text.as_bytes())));
    check(|text| Rc::<[u8]>::from(text.as_bytes()));
    check(|text| Arc::<[u8]>::from(text.as_bytes()));
    check(|text| Box::new(Vec::from(text.as_bytes())));
    check(|text| Rc::new(Vec::from(text.as_bytes())));
    check(|text| Arc::new(Vec::from(text.as_bytes())));

    check(|text| text);
    check(String::from);
    check(Box::<str>::from);
    check(Cow::Borrowed);
    check(|text| Cow::<'static, str>::Owned(String::from(text)));
    check(Rc::<str>::from);
    check(Arc::<str>::from);
    check(|text| Box::new(String::from(text)));
    check(|text| Rc::new(String::from(text)));
    check(|text| Arc::new(String::from(text)));
}

#[test]
fn stringbytes_components_support_nested_wrappers_and_trait_objects() {
    // StringBytes alone must enable construction, updates, and encoding, even
    // when the wrapper has no AsRef<[u8]> implementation.
    fn check<T: StringBytes>(wrap: impl Fn(&'static str) -> T) {
        let input = FastqParts::new(
            wrap("read1"),
            Some(wrap("note")),
            wrap("AC"),
            wrap(""),
            wrap("12"),
        );
        let pointer = input.sequence.as_bytes().as_ptr();
        let mut fields = FastqField::try_from(input).unwrap();
        let borrowed = FastqField::<&[u8]>::from(&fields);
        assert_eq!(borrowed.sequence().as_ptr(), pointer);

        // Failed updates must leave both the validated bytes and their buffer
        // in place, including when the original buffer is shared.
        assert!(matches!(
            fields.use_sequence(wrap("A\n")),
            Err(FastqFieldError::LineFeed { field: "sequence" })
        ));
        assert!(matches!(
            fields.use_sequence_and_quality(wrap("A"), wrap("12")),
            Err(FastqFieldError::LengthMismatch(_))
        ));
        assert_eq!(fields.sequence().as_bytes(), b"AC");
        assert_eq!(fields.sequence().as_bytes().as_ptr(), pointer);
        fields.use_id(wrap("read2")).unwrap();
        fields.use_description(Some(wrap("updated"))).unwrap();
        fields.use_separator(wrap("read2")).unwrap();
        fields.use_sequence(wrap("TG")).unwrap();
        fields.use_quality(wrap("34")).unwrap();
        fields
            .use_sequence_and_quality(wrap("T"), wrap("3"))
            .unwrap();

        let expected = b"@read2 updated\nT\n+read2\n3\n";
        let mut encoded = Vec::new();
        fields.encode_into_vec(&mut encoded);
        assert_eq!(encoded, expected);
        let mut written = Vec::new();
        fields.write(&mut written).unwrap();
        assert_eq!(written, expected);

        let record = FastqRecord::<Vec<u8>>::from(fields.as_ref());
        assert_eq!(
            record.parts(),
            FastqField::<&[u8]>::from(&fields).into_parts()
        );
        encoded.clear();
        record.encode_into_buf(&mut encoded);
        assert_eq!(encoded, expected);

        // The same generic constructor still rejects invalid original parts.
        let invalid = FastqParts::new(wrap("read 1"), None, wrap("AC"), wrap(""), wrap("12"));
        assert!(matches!(
            FastqField::try_from(invalid),
            Err(FastqFieldError::InvalidId)
        ));
    }

    check(|text| Rc::new(String::from(text)));
    check(|text| Arc::new(Box::new(String::from(text))));
    check(|text| Rc::new(Arc::new(Box::new(text.as_bytes().to_vec()))));
    check(|text| Box::new(String::from(text)) as Box<dyn StringBytes>);
    check(|text| Rc::new(String::from(text)) as Rc<dyn StringBytes>);
    check(|text| Arc::new(text.as_bytes().to_vec()) as Arc<dyn StringBytes>);
}

#[test]
fn fields_can_borrow_mutable_byte_components_without_copying() {
    let mut id = *b"read1";
    let mut sequence = *b"AC";
    let mut separator = [];
    let mut quality = *b"12";
    let pointer = sequence.as_ptr();
    let fields = FastqField::try_from(FastqParts::new(
        id.as_mut_slice(),
        None,
        sequence.as_mut_slice(),
        separator.as_mut_slice(),
        quality.as_mut_slice(),
    ))
    .unwrap();
    assert_eq!(fields.sequence().as_ptr(), pointer);
    let mut encoded = Vec::new();
    fields.encode_into_vec(&mut encoded);
    assert_eq!(encoded, b"@read1\nAC\n+\n12\n");
    drop(fields);
    assert_eq!(sequence, *b"AC");
}

#[test]
fn parts_conversion_accepts_arrays_and_checks_their_contents() {
    fn check<const N: usize>() {
        let parts = FastqParts::new([b'r'; N], Some([b'd'; N]), [b'A'; N], [b'r'; N], [b'!'; N]);
        let borrowed: FastqField<&[u8; N]> = FastqField::try_from(parts.as_ref()).unwrap();
        assert!(std::ptr::eq(*borrowed.sequence(), &parts.sequence));
        assert_eq!(borrowed.into_parts(), parts.as_ref());
        let fields: FastqField<[u8; N]> = FastqField::try_from(parts.clone()).unwrap();
        assert_eq!(fields.into_parts(), parts);

        if N > 0 {
            let mut invalid = parts.clone();
            invalid.sequence[0] = b'\n';
            assert!(matches!(
                FastqField::try_from(invalid.as_ref()),
                Err(FastqFieldError::LineFeed { field: "sequence" })
            ));
            assert!(matches!(
                FastqField::try_from(invalid),
                Err(FastqFieldError::LineFeed { field: "sequence" })
            ));
            let mut invalid = parts;
            invalid.id[0] = b' ';
            assert!(matches!(
                FastqField::try_from(invalid),
                Err(FastqFieldError::InvalidId)
            ));
        }
    }

    check::<0>();
    check::<1>();
    check::<2>();
    check::<64>();
}

#[test]
fn record_parts_conversion_accepts_owned_text_for_each_byte_storage() {
    fn check<U: StringBytes + AsRef<[u8]>>()
    where
        FastqRecord<U>: From<FastqField<String>>,
    {
        let record = FastqRecord::<U>::try_from(FastqParts {
            id: String::from("read1"),
            description: Some(String::from("description")),
            sequence: String::from("AC"),
            separator: String::from("read1"),
            quality: String::from("12"),
        })
        .unwrap();
        assert_eq!(record.parts(), parts());
        assert_eq!(
            FastqRecord::<&[u8]>::from(&record).into_parts().0,
            b"read1descriptionACread112"
        );
        let (_, ranges) = record.into_parts();
        assert_eq!(
            ranges,
            FastqParts {
                id: 0..5,
                description: Some(5..16),
                sequence: 16..18,
                separator: 18..23,
                quality: 23..25,
            }
        );
    }

    check::<Vec<u8>>();
    check::<Box<[u8]>>();
    check::<Bytes>();
    check::<BytesMut>();
}

#[test]
fn parts_conversion_validates_each_component_before_creating_fields() {
    for name in [
        "identifier",
        "description",
        "sequence",
        "separator",
        "quality",
    ] {
        let mut input = parts();
        match name {
            "identifier" => input.id = b"read1\nextra",
            "description" => input.description = Some(b"description\nextra"),
            "sequence" => input.sequence = b"A\n",
            "separator" => input.separator = b"+read1\nextra",
            "quality" => input.quality = b"1\n",
            _ => unreachable!(),
        }
        assert!(matches!(
            FastqField::try_from(input),
            Err(FastqFieldError::LineFeed { field }) if field == name
        ));
    }
    for id in [b"read1 extra".as_slice(), b"read1\textra"] {
        let mut input = parts();
        input.id = id;
        assert!(matches!(
            FastqField::try_from(input),
            Err(FastqFieldError::InvalidId)
        ));
    }
    let mut input = parts();
    input.quality = b"1";
    assert!(matches!(
        FastqField::try_from(input),
        Err(FastqFieldError::LengthMismatch(LengthMismatch {
            seq: 2,
            qual: 1
        }))
    ));
}

#[test]
fn use_methods_reject_invalid_replacements_without_changing_any_field() {
    let original = FastqRecord::<Bytes>::try_from(parts()).unwrap();
    let mut fields = original.fields();
    let before = fields.clone();
    for id in [b"read1 extra".as_slice(), b"read1\textra", b"read1\nextra"] {
        assert!(fields.use_id(id).is_err());
        assert_eq!(fields, before);
    }
    assert!(fields.use_description(Some(b"description\nextra")).is_err());
    assert_eq!(fields, before);
    for separator in [b"\n".as_slice(), b"read1\nextra", b"+\nextra"] {
        assert!(fields.use_separator(separator).is_err());
        assert_eq!(fields, before);
    }
    for (sequence, quality) in [
        (b"A".as_slice(), b"12".as_slice()),
        (b"A\n".as_slice(), b"12".as_slice()),
        (b"AC".as_slice(), b"1\n".as_slice()),
    ] {
        assert!(fields.use_sequence_and_quality(sequence, quality).is_err());
        assert_eq!(fields, before);
    }
    for sequence in [b"".as_slice(), b"A", b"ACG", b"A\n"] {
        assert!(fields.use_sequence(sequence).is_err());
        assert_eq!(fields, before);
    }
    for quality in [b"".as_slice(), b"1", b"123", b"1\n"] {
        assert!(fields.use_quality(quality).is_err());
        assert_eq!(fields, before);
    }

    fields.use_sequence(b"GT").unwrap();
    assert_eq!(*fields.sequence(), b"GT");
    assert_eq!(fields.quality(), before.quality());
    fields.use_quality(b"34").unwrap();
    assert_eq!(*fields.quality(), b"34");
    assert_eq!(*fields.sequence(), b"GT");

    fields.use_id(b"updated").unwrap();
    fields.use_description(Some(b"changed")).unwrap();
    fields.use_separator(b"updated").unwrap();
    fields.use_sequence_and_quality(b"C", b"2").unwrap();
    let updated = FastqRecord::<Bytes>::from(fields);
    assert_eq!(
        FastqRecord::<&[u8]>::from(&updated).into_parts().0,
        b"updatedchangedCupdated2"
    );
    assert_eq!(original.parts(), parts());
}

#[test]
fn use_methods_accept_owned_components_and_equal_length_empty_sequence_quality() {
    let mut fields = FastqField::try_from(FastqParts {
        id: b"read1".to_vec(),
        description: Some(b"description".to_vec()),
        sequence: b"AC".to_vec(),
        separator: b"read1".to_vec(),
        quality: b"12".to_vec(),
    })
    .unwrap();
    fields.use_id(Vec::new()).unwrap();
    fields.use_description(None).unwrap();
    fields.use_separator(b"".to_vec()).unwrap();
    let sequence = "é".as_bytes().to_vec();
    let sequence_pointer = sequence.as_ptr();
    fields.use_sequence(sequence).unwrap();
    assert_eq!(fields.sequence().as_ptr(), sequence_pointer);
    assert_eq!(fields.quality(), b"12");
    let quality = b"34".to_vec();
    let quality_pointer = quality.as_ptr();
    fields.use_quality(quality).unwrap();
    assert_eq!(fields.quality().as_ptr(), quality_pointer);
    assert_eq!(fields.sequence(), "é".as_bytes());
    fields
        .use_sequence_and_quality(Vec::new(), Vec::new())
        .unwrap();
    fields.use_sequence(Vec::new()).unwrap();
    fields.use_quality(Vec::new()).unwrap();
    assert_eq!(
        FastqRecord::<&[u8]>::from(&FastqRecord::<Bytes>::from(fields))
            .into_parts()
            .0,
        b""
    );
}

#[test]
fn editing_parts_requires_revalidation_and_leaves_validated_fields_unchanged() {
    let fields = FastqField::try_from(FastqParts {
        id: b"read1".to_vec(),
        description: Some(b"description".to_vec()),
        sequence: b"AC".to_vec(),
        separator: b"read1".to_vec(),
        quality: b"12".to_vec(),
    })
    .unwrap();
    assert_eq!(fields.id(), b"read1");
    assert_eq!(fields.description().unwrap(), b"description");
    assert_eq!(fields.sequence(), b"AC");
    assert_eq!(fields.separator(), b"read1");
    assert_eq!(fields.quality(), b"12");

    let mut parts = FastqField::<&[u8]>::from(&fields).into_parts();
    parts.sequence = b"C";
    assert!(matches!(
        FastqField::try_from(parts.clone()),
        Err(FastqFieldError::LengthMismatch(LengthMismatch {
            seq: 1,
            qual: 2
        }))
    ));
    parts.quality = b"2";
    let updated = FastqRecord::<Bytes>::from(FastqField::try_from(parts).unwrap());
    assert_eq!(
        FastqRecord::<&[u8]>::from(&updated).into_parts().0,
        b"read1descriptionCread12"
    );
    assert_eq!(updated.fields().into_parts(), updated.parts());
    assert_eq!(fields.sequence(), b"AC");
    assert_eq!(fields.quality(), b"12");

    let mut owned_parts = fields.into_parts();
    owned_parts.separator.push(b'\n');
    assert!(matches!(
        FastqField::try_from(owned_parts),
        Err(FastqFieldError::LineFeed { field: "separator" })
    ));
}

#[test]
fn separator_content_is_retained_without_trimming_or_inserting_a_marker() {
    for separator in [b"".as_slice(), b"read1", b"+read1", b"++", b" \tread1 \t"] {
        let input = FastqParts::new(b"read1".as_slice(), None, b"AC", separator, b"12");
        assert_eq!(input.separator, separator);
        let fields = FastqField::try_from(input.clone()).unwrap();
        assert_eq!(*fields.separator(), separator);
        assert_eq!(fields.separator().as_ptr(), separator.as_ptr());
        assert_eq!(fields.into_parts(), input);

        let mut fields = FastqField::try_from(parts()).unwrap();
        fields.use_separator(separator).unwrap();
        assert_eq!(*fields.separator(), separator);
        assert_eq!(fields.separator().as_ptr(), separator.as_ptr());
        let record = FastqRecord::<Vec<u8>>::from(fields);
        let expected = [
            b"@read1 description\nAC\n+".as_slice(),
            separator,
            b"\n12\n",
        ]
        .concat();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&record).into_parts().0,
            [b"read1descriptionAC".as_slice(), separator, b"12"].concat()
        );
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
