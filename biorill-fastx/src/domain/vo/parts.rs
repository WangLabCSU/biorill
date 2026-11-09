use std::ops::Deref;

/// FASTQ fields that can be edited directly.
///
/// Creating or editing parts does not check their contents. Convert them to
/// [`crate::FastqField`] to validate the field format and sequence/quality lengths.
/// Use [`crate::FastqField::into_parts`] to take fields out for further editing.
///
/// ```
/// use biorill_fastx::{FastqField, FastqParts};
///
/// let mut parts = FastqParts::new("read1", None, "AC", "", "12");
/// parts.description = Some("UMI:AC");
/// parts.sequence = "C";
/// parts.quality = "2";
/// let fields = FastqField::try_from(parts)?;
/// assert_eq!(*fields.sequence(), "C");
/// # Ok::<(), biorill_fastx::FastqFieldError>(())
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqParts<T> {
    /// Identifier, excluding the leading `@` and description.
    pub id: T,
    /// Description, excluding the whitespace that separates it from the identifier.
    pub description: Option<T>,
    /// Sequence.
    pub sequence: T,
    /// Content after the separator line's `+`. Use an empty value for a bare `+` line.
    pub separator: T,
    /// Quality scores corresponding to the sequence.
    pub quality: T,
}

impl<T> FastqParts<T> {
    /// Creates parts from the supplied values without validation.
    /// Supply `separator` without the line's leading `+`; values are kept as-is.
    pub fn new(id: T, description: Option<T>, sequence: T, separator: T, quality: T) -> Self {
        Self {
            id,
            description,
            sequence,
            separator,
            quality,
        }
    }

    /// Borrows each field, returning `FastqParts<&T>`.
    ///
    /// ```
    /// use biorill_fastx::FastqParts;
    ///
    /// let parts = FastqParts::new(
    ///     String::from("read1"), None, String::from("AC"),
    ///     String::new(), String::from("12"),
    /// );
    /// let borrowed: FastqParts<&String> = parts.as_ref();
    /// assert_eq!(borrowed.id, "read1");
    /// ```
    pub fn as_ref(&self) -> FastqParts<&T> {
        FastqParts {
            id: &self.id,
            description: self.description.as_ref(),
            sequence: &self.sequence,
            separator: &self.separator,
            quality: &self.quality,
        }
    }

    /// Borrows each field through [`Deref`].
    /// For example, `FastqParts<String>` becomes `FastqParts<&str>`.
    ///
    /// ```
    /// use biorill_fastx::FastqParts;
    ///
    /// let parts = FastqParts::new(
    ///     String::from("read1"), None, String::from("AC"),
    ///     String::new(), String::from("12"),
    /// );
    /// let text: FastqParts<&str> = parts.as_deref();
    /// assert_eq!(text.sequence, "AC");
    /// ```
    pub fn as_deref(&self) -> FastqParts<&T::Target>
    where
        T: Deref,
    {
        FastqParts {
            id: self.id.deref(),
            description: self.description.as_deref(),
            sequence: self.sequence.deref(),
            separator: self.separator.deref(),
            quality: self.quality.deref(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::ops::Deref;

    #[cfg(feature = "bytes")]
    use bytes::{Bytes, BytesMut};

    use super::FastqParts;

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
        #[cfg(feature = "bytes")]
        check(Bytes::from);
        #[cfg(feature = "bytes")]
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
}
