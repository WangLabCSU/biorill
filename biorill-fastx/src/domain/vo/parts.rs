use std::ops::Deref;

/// Editable FASTQ components, without validation guarantees.
///
/// Convert into [`crate::FastqField`] to check their format and sequence/quality
/// lengths. Use [`crate::FastqField::into_parts`] to recover editable components.
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
    /// Identifier without the header's leading `@` or separating space/tab.
    pub id: T,
    /// Optional description without the identifier or separating whitespace.
    pub description: Option<T>,
    /// Sequence bytes in input order.
    pub sequence: T,
    /// Content after the separator line's leading `+`; empty for a bare `+` line.
    /// Any `+` in this value is part of the content, not a format marker.
    pub separator: T,
    /// Quality bytes in the same order as the sequence.
    pub quality: T,
}

impl<T> FastqParts<T> {
    /// Create editable components without checking their format or lengths.
    /// `separator` is interpreted as content after the line's `+`; no prefix is removed.
    pub fn new(id: T, description: Option<T>, sequence: T, separator: T, quality: T) -> Self {
        Self {
            id,
            description,
            sequence,
            separator,
            quality,
        }
    }

    /// Borrow each component, converting `&FastqParts<T>` to `FastqParts<&T>`.
    /// Does not copy or validate the components.
    ///
    /// ```
    /// use biorill_fastx::FastqParts;
    ///
    /// let parts = FastqParts::new(1, None, 2, 3, 4);
    /// let borrowed: FastqParts<&i32> = parts.as_ref();
    /// assert!(std::ptr::eq(borrowed.id, &parts.id));
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

    /// Borrow each component's dereference target without validating it.
    /// Dereferences one level: `String` becomes `&str`, `Vec<u8>` becomes `&[u8]`,
    /// and `&String` becomes `&String`.
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
    /// assert_eq!(text.sequence.as_ptr(), parts.sequence.as_ptr());
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
}
