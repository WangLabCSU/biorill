use std::fmt::{self, Write};

use bytes::Bytes;

use biorill_fastx::{FastqField, FastqFieldError, FastqParts, FastqRecord};

#[test]
fn record_text_views_validate_selected_fields_and_reject_split_characters() {
    let record = FastqRecord::try_from((
        b"\xffread1AC12\xff".as_slice(),
        FastqParts::new(1..6, None, 6..8, 8..8, 8..10),
    ))
    .unwrap();
    for fields in [record.fields(), FastqField::from(&record)] {
        let text: FastqField<&str> = fields.try_into().unwrap();
        assert_eq!(text.to_string(), "@read1\nAC\n+\n12\n");
        assert_eq!(text.id().as_ptr(), record.id().as_ptr());
        assert_eq!(text.sequence().as_ptr(), record.sequence().as_ptr());
        assert_eq!(text.quality().as_ptr(), record.quality().as_ptr());
    }

    let record = FastqRecord::try_from((
        String::from("é"),
        FastqParts::new(0..0, None, 0..1, 0..0, 1..2),
    ))
    .unwrap();
    let text: FastqRecord<&str> = (&record).try_into().unwrap();
    let fields = text.fields();
    assert_eq!(*fields.sequence(), &[0xc3]);
    assert_eq!(*fields.quality(), &[0xa9]);
    assert!(matches!(
        FastqField::<&str>::try_from(fields),
        Err(FastqFieldError::InvalidUtf8Slice {
            field: "sequence",
            ..
        })
    ));
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
            let owned = FastqField::<String>::from(borrowed.clone());
            let record = FastqRecord::<Vec<u8>>::from(borrowed.clone());
            let mut expected = Vec::new();
            record.write(&mut expected).unwrap();

            assert_eq!(borrowed.to_string().as_bytes(), expected);
            assert_eq!(owned.to_string().as_bytes(), expected);
            assert_eq!(owned.as_deref().to_string().as_bytes(), expected);

            let bytes = FastqField::<Bytes>::from(owned);
            let checked: FastqField<&str> = (&bytes).try_into().unwrap();
            assert_eq!(checked.to_string().as_bytes(), expected);
            assert_eq!(checked.sequence().as_ptr(), bytes.sequence().as_ptr());
        }
    }
}

#[test]
fn custom_parts_display_through_validated_text_views() {
    struct Text(&'static str);

    impl AsRef<str> for Text {
        fn as_ref(&self) -> &str {
            self.0
        }
    }

    impl AsRef<[u8]> for Text {
        fn as_ref(&self) -> &[u8] {
            self.0.as_bytes()
        }
    }

    for description in [None, Some(""), Some("note")] {
        let parts = FastqParts {
            id: Text("read1"),
            description: description.map(Text),
            sequence: Text("AC"),
            separator: Text(""),
            quality: Text("12"),
        };
        let fields = FastqField::try_from(FastqParts::<&str>::new(
            parts.id.as_ref(),
            parts.description.as_ref().map(AsRef::as_ref),
            parts.sequence.as_ref(),
            parts.separator.as_ref(),
            parts.quality.as_ref(),
        ))
        .unwrap();
        let header = description.map_or(String::from("@read1"), |text| format!("@read1 {text}"));
        assert_eq!(
            fields.as_deref().to_string(),
            format!("{header}\nAC\n+\n12\n")
        );
    }
}

#[test]
fn display_requires_only_display_on_the_components() {
    // This helper must compile without Debug or AsRef bounds on T.
    fn format<T: fmt::Display>(fields: &FastqField<T>) -> (String, String) {
        (fields.to_string(), format!("{fields:#}"))
    }

    for (description, header) in [(None, "@r"), (Some(""), "@r "), (Some("note"), "@r note")] {
        let fields =
            FastqField::try_from(FastqParts::new("r", description, "AC", "", "12")).unwrap();
        let expected = format!("{header}\nAC\n+\n12\n");
        assert_eq!(format(&fields), (expected.clone(), expected.clone()));
        let owned = FastqField::<String>::from(fields);
        assert_eq!(format(&owned), (expected.clone(), expected));
    }
}

#[test]
fn debug_requires_only_debug_on_the_components() {
    // This helper must compile without AsRef bounds on T.
    fn format<T: fmt::Debug>(fields: &FastqField<T>) -> (String, String) {
        (format!("{fields:?}"), format!("{fields:#?}"))
    }

    let fields = FastqField::try_from(FastqParts::new("r", None, "AC", "", "12")).unwrap();
    let (compact, pretty) = format(&fields);
    assert_eq!(
        compact,
        "FastqField { id: \"r\", description: None, sequence: \"AC\", separator: \"\", quality: \"12\" }",
    );
    assert!(pretty.starts_with("FastqField {\n"));
    assert!(pretty.contains("    id: \"r\",\n"));
    assert!(pretty.contains("    description: None,\n"));
}

#[test]
fn display_propagates_output_errors_and_stops_writing() {
    struct LimitedWriter {
        remaining: usize,
        failed: bool,
    }

    impl Write for LimitedWriter {
        fn write_str(&mut self, text: &str) -> fmt::Result {
            assert!(!self.failed, "formatting continued after a write error");
            if text.len() > self.remaining {
                self.failed = true;
                return Err(fmt::Error);
            }
            self.remaining -= text.len();
            Ok(())
        }
    }

    let fields = FastqField::try_from(FastqParts {
        id: "read1",
        description: Some("说明"),
        sequence: "AC",
        separator: "read1",
        quality: "12",
    })
    .unwrap();
    let owned = FastqField::<String>::from(fields.clone());
    for value in [&fields as &dyn fmt::Display, &owned] {
        let expected = value.to_string();
        for remaining in 0..=expected.len() {
            let mut writer = LimitedWriter {
                remaining,
                failed: false,
            };
            let result = write!(writer, "{value}");
            assert_eq!(result.is_ok(), remaining == expected.len());
            assert_eq!(writer.failed, remaining < expected.len());
        }
    }
    assert_eq!(fields.to_string(), "@read1 说明\nAC\n+read1\n12\n");
}
