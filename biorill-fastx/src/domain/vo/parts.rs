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
