use std::fmt;
use std::io::{self, Write};
use std::ops::Deref;
use std::str::Utf8Error;
use std::string::FromUtf8Error;

use bytes::BufMut;
use stringbytes::StringBytes;
use thiserror::Error;

use super::parts::FastqParts;
use super::record::FastqRecord;

/// Validated FASTQ fields with checked updates.
///
/// Construct from [`FastqParts`] with [`Self::try_from`]. Components must implement
/// [`StringBytes`], which supports byte and text storage, references, shared
/// pointers, and nested wrappers. Validation scans the fields and moves the
/// supplied components without copying them.
/// The `set_*` methods validate new values before changing any field.
/// Sequence and quality lengths must match; alphabets and encodings are not checked.
/// The separator contains only the content after the line's leading `+`.
/// Use an empty value for a bare `+` line; supplied content is never trimmed.
///
/// [`fmt::Display`] is available when `T: fmt::Display`. Each component uses its
/// own `Display` implementation in a four-line FASTQ layout with LF endings,
/// including the final newline, a leading `+` on the separator line, and a space
/// before any description. `&str` and `String` fields produce FASTQ text.
/// To display byte fields as FASTQ text, use [`Self::try_str_view`] to check UTF-8
/// and borrow the fields as text.
///
/// [`fmt::Debug`] uses the derived field structure: `{:?}` is compact and `{:#?}`
/// is pretty-printed.
///
/// # Examples
///
/// ```
/// use biorill_fastx::{FastqField, FastqParts};
///
/// let mut fields = FastqField::try_from(FastqParts::new(
///     "read1",
///     None,
///     "AC",
///     "",
///     "12",
/// ))?;
/// fields.set_description(Some("UMI:AC"))?;
/// fields.set_sequence_and_quality("C", "2")?;
/// assert_eq!(fields.to_string(), "@read1 UMI:AC\nC\n+\n2\n");
/// # Ok::<(), biorill_fastx::FastqFieldError>(())
/// ```
///
/// ```compile_fail,E0616
/// use biorill_fastx::{FastqField, FastqParts};
///
/// let mut fields = FastqField::try_from(FastqParts::new(
///     "read1",
///     None,
///     "AC",
///     "",
///     "12",
/// )).unwrap();
/// fields.sequence = "A";
/// ```
///
/// Implementing `AsRef<[u8]>` alone does not enable construction from custom
/// components. Borrow each custom component as a byte slice and construct
/// `FastqParts<&[u8]>` for validation instead.
///
/// ```compile_fail,E0277
/// use biorill_fastx::{FastqField, FastqParts};
///
/// struct Custom(&'static [u8]);
/// impl AsRef<[u8]> for Custom {
///     fn as_ref(&self) -> &[u8] { self.0 }
/// }
///
/// let parts = FastqParts::new(
///     Custom(b"read1"), None, Custom(b"AC"), Custom(b""), Custom(b"12"),
/// );
/// let fields = FastqField::<Custom>::try_from(parts);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqField<T> {
    id: T,
    description: Option<T>,
    sequence: T,
    separator: T,
    quality: T,
}

impl<T> FastqField<T> {
    /// The identifier without the leading `@` or description.
    pub fn id(&self) -> &T {
        &self.id
    }

    /// The optional description without the identifier or separating whitespace.
    pub fn description(&self) -> Option<&T> {
        self.description.as_ref()
    }

    /// The sequence in input order.
    pub fn sequence(&self) -> &T {
        &self.sequence
    }

    /// The content after the separator line's leading `+`, possibly empty.
    pub fn separator(&self) -> &T {
        &self.separator
    }

    /// The quality scores in the same order as the sequence.
    pub fn quality(&self) -> &T {
        &self.quality
    }

    /// Map the components, stopping at the first error, without validating the result.
    ///
    /// # Safety
    /// Every successful call to `map` must preserve the input component's
    /// validated byte sequence.
    unsafe fn try_map_fields<U, E>(
        self,
        map: impl Fn(&'static str, T) -> Result<U, E>,
    ) -> Result<FastqField<U>, E> {
        Ok(FastqField {
            id: map("identifier", self.id)?,
            description: self
                .description
                .map(|value| map("description", value))
                .transpose()?,
            sequence: map("sequence", self.sequence)?,
            separator: map("separator", self.separator)?,
            quality: map("quality", self.quality)?,
        })
    }

    /// Map the components without validating the result.
    ///
    /// # Safety
    /// Every call to `map` must preserve the input component's validated byte sequence.
    unsafe fn map_fields<U>(self, map: impl Fn(T) -> U) -> FastqField<U> {
        FastqField {
            id: map(self.id),
            description: self.description.map(&map),
            sequence: map(self.sequence),
            separator: map(self.separator),
            quality: map(self.quality),
        }
    }

    /// Borrow each component, converting `&FastqField<T>` to `FastqField<&T>`.
    /// Does not copy the components or revalidate them.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     String::from("read1"),
    ///     Some(String::from("note")),
    ///     String::from("AC"),
    ///     String::from(""),
    ///     String::from("12"),
    /// ))?;
    /// let borrowed: FastqField<&String> = fields.as_ref();
    /// assert!(std::ptr::eq(*borrowed.id(), fields.id()));
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn as_ref(&self) -> FastqField<&T> {
        FastqField {
            id: &self.id,
            description: self.description.as_ref(),
            sequence: &self.sequence,
            separator: &self.separator,
            quality: &self.quality,
        }
    }

    /// Consume the validated fields and return their editable components.
    pub fn into_parts(self) -> FastqParts<T> {
        FastqParts {
            id: self.id,
            description: self.description,
            sequence: self.sequence,
            separator: self.separator,
            quality: self.quality,
        }
    }
}

impl<T: StringBytes> FastqField<T> {
    /// Borrow each component's dereference target without copying or revalidating it.
    /// Dereferences one level: `String` becomes `&str`, `Vec<u8>` becomes `&[u8]`,
    /// `Rc<String>` becomes `&String`, and `&String` stays `&String`.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     b"read1".to_vec(), None, b"AC".to_vec(), Vec::new(), b"12".to_vec(),
    /// ))?;
    /// let slices: FastqField<&[u8]> = fields.as_deref();
    /// assert_eq!(*slices.sequence(), b"AC");
    /// assert_eq!(slices.sequence().as_ptr(), fields.sequence().as_ptr());
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn as_deref(&self) -> FastqField<&T::Target>
    where
        T: Deref,
    {
        // SAFETY: StringBytes guarantees that Deref::Target implements StringBytes
        // and exposes the same byte slice, preserving each component's validated bytes.
        unsafe { self.as_ref().map_fields(Deref::deref) }
    }

    /// Return a view of these fields with each component borrowed as a byte slice.
    /// Does not allocate, copy, or scan the field contents.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     "read1", None, "AC", "", "12",
    /// ))?;
    /// let view = fields.byte_view();
    /// assert_eq!(*view.sequence(), b"AC");
    /// assert_eq!(view.sequence().as_ptr(), fields.sequence().as_ptr());
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn byte_view(&self) -> FastqField<&[u8]> {
        // SAFETY: StringBytes::as_bytes exposes each component's validated bytes.
        unsafe { self.as_ref().map_fields(T::as_bytes) }
    }

    /// Check UTF-8 where needed and borrow each component as text without copying.
    /// Components that already guarantee UTF-8 skip validation.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::InvalidUtf8Slice`] for the first invalid field.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     b"r1".to_vec(), None, b"AC".to_vec(), Vec::new(), b"12".to_vec(),
    /// ))?;
    /// let text = fields.try_str_view()?;
    /// assert_eq!(*text.sequence(), "AC");
    /// assert_eq!(text.sequence().as_ptr(), fields.sequence().as_ptr());
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn try_str_view(&self) -> Result<FastqField<&str>, FastqFieldError> {
        // SAFETY: Successful StringBytes text views preserve the input bytes.
        unsafe {
            self.as_ref().try_map_fields(|field, value| {
                T::try_borrow_str(value)
                    .map_err(|source| FastqFieldError::InvalidUtf8Slice { field, source })
            })
        }
    }

    /// Set `id` as the identifier, without its leading `@` or description.
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing the fields if `id` contains
    /// LF, a space, or a tab.
    pub fn set_id(&mut self, id: T) -> Result<(), FastqFieldError> {
        validate_line("identifier", id.as_bytes())?;
        validate_id_delimiters(id.as_bytes())?;
        self.id = id;
        Ok(())
    }

    /// Set the description, or remove it with `None`.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::LineFeed`] without changing the fields if the
    /// description contains LF.
    pub fn set_description(&mut self, description: Option<T>) -> Result<(), FastqFieldError> {
        if let Some(description) = &description {
            validate_line("description", description.as_bytes())?;
        }
        self.description = description;
        Ok(())
    }

    /// Set the content after the separator line's leading `+`.
    /// An empty value produces a bare `+` line; supplied content is never trimmed.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::LineFeed`] without changing the fields if the
    /// separator content contains LF.
    pub fn set_separator(&mut self, separator: T) -> Result<(), FastqFieldError> {
        validate_line("separator", separator.as_bytes())?;
        self.separator = separator;
        Ok(())
    }

    /// Set a new sequence with the same byte length as the current quality scores.
    /// To change both lengths, use [`Self::set_sequence_and_quality`].
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing the fields if the sequence
    /// contains LF or its byte length differs from the current quality length.
    pub fn set_sequence(&mut self, sequence: T) -> Result<(), FastqFieldError> {
        validate_line("sequence", sequence.as_bytes())?;
        validate_lengths(sequence.as_bytes(), self.quality.as_bytes())?;
        self.sequence = sequence;
        Ok(())
    }

    /// Set new quality scores with the same byte length as the current sequence.
    /// To change both lengths, use [`Self::set_sequence_and_quality`].
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing the fields if the quality
    /// contains LF or its byte length differs from the current sequence length.
    pub fn set_quality(&mut self, quality: T) -> Result<(), FastqFieldError> {
        validate_line("quality", quality.as_bytes())?;
        validate_lengths(self.sequence.as_bytes(), quality.as_bytes())?;
        self.quality = quality;
        Ok(())
    }

    /// Set the sequence and quality together. Empty values are accepted.
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing either field if a replacement
    /// contains LF or the replacement lengths differ.
    pub fn set_sequence_and_quality(
        &mut self,
        sequence: T,
        quality: T,
    ) -> Result<(), FastqFieldError> {
        validate_line("sequence", sequence.as_bytes())?;
        validate_line("quality", quality.as_bytes())?;
        validate_lengths(sequence.as_bytes(), quality.as_bytes())?;
        self.sequence = sequence;
        self.quality = quality;
        Ok(())
    }

    /// Write four FASTQ lines with LF endings and a space before any description.
    /// Field bytes are written without checking their encoding.
    ///
    /// # Errors
    /// Returns the writer's error if the fields cannot be written in full.
    /// A prefix may already have been written.
    pub fn write<W: Write + ?Sized>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(b"@")?;
        writer.write_all(self.id.as_bytes())?;
        if let Some(description) = &self.description {
            writer.write_all(b" ")?;
            writer.write_all(description.as_bytes())?;
        }
        writer.write_all(b"\n")?;
        writer.write_all(self.sequence.as_bytes())?;
        writer.write_all(b"\n+")?;
        writer.write_all(self.separator.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(self.quality.as_bytes())?;
        writer.write_all(b"\n")
    }

    /// Encode the fields and append them to `buf`, reserving space for the
    /// complete record first. Uses LF endings and a space before any description.
    /// Copies all field bytes; growing `buf` may allocate.
    pub fn encode_into_vec(&self, buf: &mut Vec<u8>) {
        let amt = self.id.as_bytes().len()
            + self
                .description
                .as_ref()
                .map_or(0, |desc| desc.as_bytes().len() + 1)
            + self.sequence.as_bytes().len()
            + self.separator.as_bytes().len()
            + self.quality.as_bytes().len()
            + 6;
        buf.reserve(amt);
        self.encode_into_buf(buf);
    }

    /// Encode four FASTQ lines and append them to `buf`, using LF endings and a
    /// space before any description. Field bytes are copied into `buf`.
    ///
    /// # Panics
    /// Panics if `buf` cannot accept the fields. A prefix may already have been appended.
    pub fn encode_into_buf<B: BufMut + ?Sized>(&self, buf: &mut B) {
        buf.put_slice(b"@");
        buf.put_slice(self.id.as_bytes());
        if let Some(description) = &self.description {
            buf.put_slice(b" ");
            buf.put_slice(description.as_bytes());
        }
        buf.put_slice(b"\n");
        buf.put_slice(self.sequence.as_bytes());
        buf.put_slice(b"\n+");
        buf.put_slice(self.separator.as_bytes());
        buf.put_slice(b"\n");
        buf.put_slice(self.quality.as_bytes());
        buf.put_slice(b"\n");
    }
}

impl FastqField<Vec<u8>> {
    /// Consume the fields and convert each component into a [`Vec<u8>`].
    /// Reuses storage where possible and copies otherwise, following
    /// [`StringBytes::into_bytes`]. Does not validate UTF-8.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     String::from("r1"), None, String::from("AC"),
    ///     String::new(), String::from("12"),
    /// ))?;
    /// let pointer = fields.sequence().as_ptr();
    /// let bytes = FastqField::<Vec<u8>>::from_fields(fields);
    /// assert_eq!(bytes.sequence(), b"AC");
    /// assert_eq!(bytes.sequence().as_ptr(), pointer);
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn from_fields<T: StringBytes>(fields: FastqField<T>) -> Self {
        // SAFETY: StringBytes::into_bytes preserves each component's bytes.
        unsafe { fields.map_fields(StringBytes::into_bytes) }
    }
}

impl FastqField<String> {
    /// Consume the fields and convert each component into a [`String`].
    /// Reuses storage where possible and copies otherwise, following
    /// [`StringBytes::try_into_string`]. Checks UTF-8 where needed.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::InvalidUtf8`] for the first invalid field,
    /// retaining its name and rejected bytes. The other components are dropped.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     b"r1".to_vec(), None, b"AC".to_vec(), Vec::new(), b"12".to_vec(),
    /// ))?;
    /// let copied = FastqField::<String>::try_from_fields(fields.as_ref())?;
    /// let owned = FastqField::<String>::try_from_fields(fields)?;
    /// assert_eq!(copied, owned);
    /// assert_eq!(owned.sequence(), "AC");
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn try_from_fields<T: StringBytes>(fields: FastqField<T>) -> Result<Self, FastqFieldError> {
        // SAFETY: Successful StringBytes::try_into_string conversions preserve the input bytes.
        unsafe {
            fields.try_map_fields(|field, value| {
                value
                    .try_into_string()
                    .map_err(|source| FastqFieldError::InvalidUtf8 { field, source })
            })
        }
    }
}

impl<'a, T: StringBytes> From<&'a FastqRecord<T>> for FastqField<&'a [u8]> {
    /// Borrow the validated fields without copying or revalidating them.
    /// Updating the returned fields does not change the source record.
    fn from(record: &'a FastqRecord<T>) -> Self {
        // The record already validates these ranges and fields. StringBytes
        // keeps the borrowed bytes unchanged through shared access.
        Self {
            id: record.id(),
            description: record.description(),
            sequence: record.sequence(),
            separator: record.separator(),
            quality: record.quality(),
        }
    }
}

impl<T: fmt::Display> fmt::Display for FastqField<T> {
    /// Format four FASTQ lines with LF endings and a space before any description.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "@{}", self.id)?;
        if let Some(description) = &self.description {
            write!(formatter, " {description}")?;
        }
        write!(
            formatter,
            "\n{}\n+{}\n{}\n",
            self.sequence, self.separator, self.quality,
        )
    }
}

/// Sequence and quality have different byte lengths.
#[derive(Debug, Error)]
#[error("sequence and quality lengths do not match ({seq} vs {qual})")]
pub struct LengthMismatch {
    /// The sequence byte length.
    pub seq: usize,
    /// The quality byte length.
    pub qual: usize,
}

/// FASTQ fields could not be validated or converted to text.
#[derive(Debug, Error)]
pub enum FastqFieldError {
    /// A field contains an LF byte that would introduce another record line.
    #[error("FASTQ {field} contains a line feed")]
    LineFeed { field: &'static str },
    /// An identifier contains a FASTQ header delimiter.
    #[error("FASTQ identifier contains a space or tab")]
    InvalidId,
    /// Sequence and quality byte lengths differ.
    #[error(transparent)]
    LengthMismatch(#[from] LengthMismatch),
    /// A field cannot be converted to UTF-8 text.
    #[error("FASTQ {field} contains invalid UTF-8: {source}")]
    InvalidUtf8 {
        /// The field containing invalid UTF-8.
        field: &'static str,
        /// The UTF-8 error and the rejected field bytes.
        #[source]
        source: FromUtf8Error,
    },
    /// A borrowed field cannot be viewed as UTF-8 text.
    #[error("FASTQ {field} contains invalid UTF-8: {source}")]
    InvalidUtf8Slice {
        /// The field containing invalid UTF-8.
        field: &'static str,
        /// The UTF-8 validation error.
        #[source]
        source: Utf8Error,
    },
}

impl<T: StringBytes> TryFrom<FastqParts<T>> for FastqField<T> {
    type Error = FastqFieldError;

    /// Validate the parts and move their components into the fields without copying.
    /// The separator is interpreted as content after the line's `+`.
    /// A leading `+` in the supplied value is retained as content.
    /// Sequence alphabets and quality encodings are not checked.
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] for a component containing LF, an identifier
    /// containing space/tab, or unequal sequence and quality lengths. Empty
    /// identifiers, separator contents and equal-length empty sequence/quality
    /// components are accepted.
    fn try_from(parts: FastqParts<T>) -> Result<Self, Self::Error> {
        // Validate borrowed bytes, then move the original components unchanged.
        // StringBytes keeps their byte views stable while the fields own them.
        let fields = FastqParts::new(
            parts.id.as_bytes(),
            parts.description.as_ref().map(StringBytes::as_bytes),
            parts.sequence.as_bytes(),
            parts.separator.as_bytes(),
            parts.quality.as_bytes(),
        );
        for (field, bytes) in [
            ("identifier", fields.id),
            ("description", fields.description.unwrap_or_default()),
            ("sequence", fields.sequence),
            ("separator", fields.separator),
            ("quality", fields.quality),
        ] {
            validate_line(field, bytes)?;
        }
        validate_id_delimiters(fields.id)?;
        validate_lengths(fields.sequence, fields.quality)?;

        Ok(Self {
            id: parts.id,
            description: parts.description,
            sequence: parts.sequence,
            separator: parts.separator,
            quality: parts.quality,
        })
    }
}

fn validate_line(field: &'static str, bytes: &[u8]) -> Result<(), FastqFieldError> {
    if bytes.contains(&b'\n') {
        return Err(FastqFieldError::LineFeed { field });
    }
    Ok(())
}

fn validate_id_delimiters(id: &[u8]) -> Result<(), FastqFieldError> {
    if id.iter().any(|byte| matches!(byte, b' ' | b'\t')) {
        return Err(FastqFieldError::InvalidId);
    }
    Ok(())
}

fn validate_lengths(sequence: &[u8], quality: &[u8]) -> Result<(), FastqFieldError> {
    if sequence.len() != quality.len() {
        return Err(LengthMismatch {
            seq: sequence.len(),
            qual: quality.len(),
        }
        .into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;
    use std::error::Error;
    use std::fmt::{self, Write};
    use std::ops::Deref;
    use std::rc::Rc;
    use std::str::Utf8Error;
    use std::string::FromUtf8Error;
    use std::sync::Arc;

    use bytes::{Bytes, BytesMut};
    use stringbytes::StringBytes;

    use super::{FastqField, FastqFieldError, LengthMismatch};
    use crate::tests::valid_parts;
    use crate::{FastqParts, FastqRecord};

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
        let borrowed = fields.byte_view().into_parts();
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
    fn byte_view_and_from_fields_preserve_binary_fields_through_nested_storage() {
        fn check<T: StringBytes>(wrap: impl Fn(&'static [u8]) -> T) {
            for description in [None, Some(b"".as_slice()), Some(b"\xff".as_slice())] {
                let expected = FastqParts {
                    description,
                    ..parts()
                };
                let source = FastqField::try_from(FastqParts::new(
                    wrap(expected.id),
                    expected.description.map(&wrap),
                    wrap(expected.sequence),
                    wrap(expected.separator),
                    wrap(expected.quality),
                ))
                .unwrap();
                let bytes = source.byte_view();
                assert_eq!(bytes.clone().into_parts(), expected);
                for (original, borrowed) in [
                    (source.id(), *bytes.id()),
                    (source.sequence(), *bytes.sequence()),
                    (source.separator(), *bytes.separator()),
                    (source.quality(), *bytes.quality()),
                ] {
                    assert_eq!(original.as_bytes().as_ptr(), borrowed.as_ptr());
                }
                assert_eq!(
                    source.description().map(|value| value.as_bytes().as_ptr()),
                    bytes.description().map(|value| value.as_ptr())
                );
                let owned = FastqField::<Vec<u8>>::from_fields(source);
                assert_eq!(owned.byte_view().into_parts(), expected);
            }
        }

        check(|bytes| bytes);
        check(Vec::from);
        check(Box::<[u8]>::from);
        check(Bytes::from_static);
        check(BytesMut::from);
        check(|bytes| Rc::new(Arc::new(Box::<[u8]>::from(bytes))));
    }

    #[test]
    fn from_fields_reuses_vector_buffers() {
        for spare in [0, 32] {
            let fields = vectors(spare);
            let pointers = check(&fields);
            let fields = FastqField::<Vec<u8>>::from_fields(fields);
            assert_eq!(check(&fields), pointers);
        }
    }

    #[test]
    fn from_fields_copies_borrowed_storage_and_preserves_optional_descriptions() {
        for description in [None, Some(b"".as_slice()), Some(b"\xff".as_slice())] {
            let fields = {
                let mut source = vectors(0);
                source.set_description(description.map(Vec::from)).unwrap();
                let fields = FastqField::<Vec<u8>>::from_fields(source.as_deref());
                assert_eq!(fields.as_deref(), source.as_deref());
                assert_ne!(fields.id().as_ptr(), source.id().as_ptr());
                fields
            };
            assert_eq!(
                fields.as_deref().into_parts(),
                FastqParts {
                    description,
                    ..parts()
                },
            );
        }
    }

    #[test]
    fn try_from_fields_retains_the_rejected_buffer() {
        let fields = vectors(32);
        let pointer = fields.sequence().as_ptr();
        let error = FastqField::<String>::try_from_fields(fields).unwrap_err();
        let FastqFieldError::InvalidUtf8 { field, source } = error else {
            panic!("expected an owned UTF-8 error");
        };
        assert_eq!(field, "sequence");
        assert_eq!(source.as_bytes(), b"A\xff");
        assert_eq!(source.as_bytes().as_ptr(), pointer);
    }

    #[test]
    fn from_fields_copies_shared_and_static_storage() {
        for wrap in [Bytes::from_static, Bytes::copy_from_slice] {
            let input = parts();
            let fields = FastqField::try_from(FastqParts::new(
                wrap(input.id),
                input.description.map(wrap),
                wrap(input.sequence),
                wrap(input.separator),
                wrap(input.quality),
            ))
            .unwrap();
            let source = check(&fields);
            let vector = FastqField::<Vec<u8>>::from_fields(fields.clone());
            assert_ne!(check(&vector).sequence, source.sequence);
            let mut owned = vector.into_parts();
            owned.sequence[0] = b'T';
            assert_eq!(fields.sequence().as_ref(), input.sequence);
        }
    }

    #[test]
    fn try_str_view_borrows_byte_storage_and_preserves_field_contents() {
        fn check<T: StringBytes + AsRef<[u8]>>(wrap: impl Fn(&str) -> T) {
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
                    let borrowed = bytes.try_str_view().unwrap();
                    let slices = bytes.byte_view();
                    let from_slices = slices.try_str_view().unwrap();
                    assert_eq!(borrowed, from_slices);
                    assert_eq!(
                        borrowed.byte_view().into_parts(),
                        slices.clone().into_parts()
                    );
                    for text in [&borrowed, &from_slices] {
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

        check(|text| text.to_owned());
        check(|text| text.as_bytes().to_vec());
        check(|text| Bytes::copy_from_slice(text.as_bytes()));
        check(|text| BytesMut::from(text));
        check(|text| text.as_bytes().to_vec().into_boxed_slice());
    }

    #[test]
    fn try_from_fields_copies_nonstatic_slices_and_outlives_the_source() {
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
                    let slices = source.byte_view();
                    let owned = FastqField::<String>::try_from_fields(slices).unwrap();
                    assert_eq!(owned.byte_view(), source.byte_view());
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
    fn try_from_fields_preserves_generic_support_and_reuses_nested_owned_storage() {
        fn check<T: StringBytes>(wrap: impl Fn(String) -> T) {
            let fields = FastqField::try_from(FastqParts::new(
                wrap(String::from("read1")),
                Some(wrap(String::from("note"))),
                wrap(String::from("AC")),
                wrap(String::from("read1")),
                wrap(String::from("12")),
            ))
            .unwrap();
            let pointers = |fields: &FastqField<T>| {
                [
                    fields.id().as_bytes().as_ptr(),
                    fields.description().unwrap().as_bytes().as_ptr(),
                    fields.sequence().as_bytes().as_ptr(),
                    fields.separator().as_bytes().as_ptr(),
                    fields.quality().as_bytes().as_ptr(),
                ]
            };
            let original = pointers(&fields);
            let copied = FastqField::<String>::try_from_fields(fields.as_ref()).unwrap();
            assert_eq!(copied.byte_view(), fields.byte_view());
            for (copied, original) in [
                copied.id(),
                copied.description().unwrap(),
                copied.sequence(),
                copied.separator(),
                copied.quality(),
            ]
            .into_iter()
            .zip(original)
            {
                assert_ne!(copied.as_ptr(), original);
            }
            let owned = FastqField::<String>::try_from_fields(fields).unwrap();
            assert_eq!(owned, copied);
            assert_eq!(
                [
                    owned.id().as_ptr(),
                    owned.description().unwrap().as_ptr(),
                    owned.sequence().as_ptr(),
                    owned.separator().as_ptr(),
                    owned.quality().as_ptr(),
                ],
                original
            );
        }

        check(std::convert::identity);
        check(String::into_bytes);
        check(Bytes::from);
        check(|text| BytesMut::from(Bytes::from(text)));
        check(|text| text.into_bytes().into_boxed_slice());
        check(|text| Rc::new(Arc::new(text.into_bytes())));
        check(|text| Arc::new(Box::new(Rc::new(text))));
    }

    #[test]
    fn text_conversions_report_each_invalid_field_and_preserve_rejected_bytes() {
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
            let fields =
                FastqField::<Vec<u8>>::from_fields(FastqField::try_from(input.clone()).unwrap());
            let error = fields.try_str_view().unwrap_err();
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
            assert_eq!(fields.byte_view().into_parts(), input);

            for result in [
                FastqField::<String>::try_from_fields(fields.as_ref()),
                FastqField::<String>::try_from_fields(fields),
            ] {
                let error = result.unwrap_err();
                let FastqFieldError::InvalidUtf8 { field, source } = &error else {
                    panic!("expected an owned UTF-8 error");
                };
                assert_eq!(*field, name);
                assert_eq!(source.as_bytes(), invalid);
                assert_eq!(source.utf8_error().valid_up_to(), valid_prefix);
                assert_eq!(source.utf8_error().error_len(), Some(1));
                assert_eq!(
                    error.source().unwrap().downcast_ref::<FromUtf8Error>(),
                    Some(source)
                );
            }
        }
    }

    #[test]
    fn parts_conversion_accepts_common_storage_and_borrows_without_copying() {
        fn check<T>(wrap: impl Fn(&'static str) -> T)
        where
            T: StringBytes + Deref,
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
                let input =
                    FastqParts::new(wrap(id), None, wrap(sequence), wrap(""), wrap(quality));
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
            let borrowed = fields.byte_view();
            assert_eq!(borrowed.sequence().as_ptr(), pointer);

            // Failed updates must leave both the validated bytes and their buffer
            // in place, including when the original buffer is shared.
            assert!(matches!(
                fields.set_sequence(wrap("A\n")),
                Err(FastqFieldError::LineFeed { field: "sequence" })
            ));
            assert!(matches!(
                fields.set_sequence_and_quality(wrap("A"), wrap("12")),
                Err(FastqFieldError::LengthMismatch(_))
            ));
            assert_eq!(fields.sequence().as_bytes(), b"AC");
            assert_eq!(fields.sequence().as_bytes().as_ptr(), pointer);
            fields.set_id(wrap("read2")).unwrap();
            fields.set_description(Some(wrap("updated"))).unwrap();
            fields.set_separator(wrap("read2")).unwrap();
            fields.set_sequence(wrap("TG")).unwrap();
            fields.set_quality(wrap("34")).unwrap();
            fields
                .set_sequence_and_quality(wrap("T"), wrap("3"))
                .unwrap();

            let expected = b"@read2 updated\nT\n+read2\n3\n";
            let mut encoded = Vec::new();
            fields.encode_into_vec(&mut encoded);
            assert_eq!(encoded, expected);
            let mut written = Vec::new();
            fields.write(&mut written).unwrap();
            assert_eq!(written, expected);

            let record = FastqRecord::<Vec<u8>>::from(fields.as_ref());
            assert_eq!(record.parts(), fields.byte_view().into_parts());
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
            let parts =
                FastqParts::new([b'r'; N], Some([b'd'; N]), [b'A'; N], [b'r'; N], [b'!'; N]);
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
    fn parts_conversion_validates_each_component_before_creating_fields() {
        for name in [
            "identifier",
            "description",
            "sequence",
            "separator",
            "quality",
        ] {
            let mut input = valid_parts();
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
            let mut input = valid_parts();
            input.id = id;
            assert!(matches!(
                FastqField::try_from(input),
                Err(FastqFieldError::InvalidId)
            ));
        }
        let mut input = valid_parts();
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
    fn set_methods_reject_invalid_replacements_without_changing_any_field() {
        let original = FastqRecord::<Bytes>::try_from(valid_parts()).unwrap();
        let mut fields = original.fields();
        let before = fields.clone();
        for id in [b"read1 extra".as_slice(), b"read1\textra", b"read1\nextra"] {
            assert!(fields.set_id(id).is_err());
            assert_eq!(fields, before);
        }
        assert!(fields.set_description(Some(b"description\nextra")).is_err());
        assert_eq!(fields, before);
        for separator in [b"\n".as_slice(), b"read1\nextra", b"+\nextra"] {
            assert!(fields.set_separator(separator).is_err());
            assert_eq!(fields, before);
        }
        for (sequence, quality) in [
            (b"A".as_slice(), b"12".as_slice()),
            (b"A\n".as_slice(), b"12".as_slice()),
            (b"AC".as_slice(), b"1\n".as_slice()),
        ] {
            assert!(fields.set_sequence_and_quality(sequence, quality).is_err());
            assert_eq!(fields, before);
        }
        for sequence in [b"".as_slice(), b"A", b"ACG", b"A\n"] {
            assert!(fields.set_sequence(sequence).is_err());
            assert_eq!(fields, before);
        }
        for quality in [b"".as_slice(), b"1", b"123", b"1\n"] {
            assert!(fields.set_quality(quality).is_err());
            assert_eq!(fields, before);
        }

        fields.set_sequence(b"GT").unwrap();
        assert_eq!(*fields.sequence(), b"GT");
        assert_eq!(fields.quality(), before.quality());
        fields.set_quality(b"34").unwrap();
        assert_eq!(*fields.quality(), b"34");
        assert_eq!(*fields.sequence(), b"GT");

        fields.set_id(b"updated").unwrap();
        fields.set_description(Some(b"changed")).unwrap();
        fields.set_separator(b"updated").unwrap();
        fields.set_sequence_and_quality(b"C", b"2").unwrap();
        let updated = FastqRecord::<Bytes>::from(fields);
        assert_eq!(
            FastqRecord::<&[u8]>::from(&updated).into_parts().0,
            b"updatedchangedCupdated2"
        );
        assert_eq!(original.parts(), valid_parts());
    }

    #[test]
    fn set_methods_accept_owned_components_and_equal_length_empty_sequence_quality() {
        let mut fields = FastqField::try_from(FastqParts {
            id: b"read1".to_vec(),
            description: Some(b"description".to_vec()),
            sequence: b"AC".to_vec(),
            separator: b"read1".to_vec(),
            quality: b"12".to_vec(),
        })
        .unwrap();
        fields.set_id(Vec::new()).unwrap();
        fields.set_description(None).unwrap();
        fields.set_separator(b"".to_vec()).unwrap();
        let sequence = "é".as_bytes().to_vec();
        let sequence_pointer = sequence.as_ptr();
        fields.set_sequence(sequence).unwrap();
        assert_eq!(fields.sequence().as_ptr(), sequence_pointer);
        assert_eq!(fields.quality(), b"12");
        let quality = b"34".to_vec();
        let quality_pointer = quality.as_ptr();
        fields.set_quality(quality).unwrap();
        assert_eq!(fields.quality().as_ptr(), quality_pointer);
        assert_eq!(fields.sequence(), "é".as_bytes());
        fields
            .set_sequence_and_quality(Vec::new(), Vec::new())
            .unwrap();
        fields.set_sequence(Vec::new()).unwrap();
        fields.set_quality(Vec::new()).unwrap();
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

        let mut parts = fields.byte_view().into_parts();
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

            let mut fields = FastqField::try_from(valid_parts()).unwrap();
            fields.set_separator(separator).unwrap();
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
            let owned = FastqField::<String>::try_from_fields(fields).unwrap();
            assert_eq!(format(&owned), (expected.clone(), expected));
        }
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
        let owned = FastqField::<String>::try_from_fields(fields.clone()).unwrap();
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
}
