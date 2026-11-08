use std::fmt;
use std::io::{self, Write};
use std::ops::Deref;
use std::str::Utf8Error;
use std::string::FromUtf8Error;

use bytes::{BufMut, Bytes, BytesMut};
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
/// The `use_*` methods validate new values before changing any field.
/// Sequence and quality lengths must match; alphabets and encodings are not checked.
/// The separator contains only the content after the line's leading `+`.
/// Use an empty value for a bare `+` line; supplied content is never trimmed.
///
/// [`fmt::Display`] is available when `T: fmt::Display`. Each component uses its
/// own `Display` implementation in a four-line FASTQ layout with LF endings,
/// including the final newline, a leading `+` on the separator line, and a space
/// before any description. `&str` and `String` fields produce FASTQ text.
/// To display byte fields as FASTQ text, first use [`TryFrom`] to borrow them as
/// `FastqField<&str>`, checking UTF-8 without copying the field bytes.
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
/// fields.use_description(Some("UMI:AC"))?;
/// fields.use_sequence_and_quality("C", "2")?;
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
///
/// Owned fields can be converted between `Vec<u8>`, [`Bytes`], [`BytesMut`] and `Box<[u8]>`
/// with [`From`]. Converting from `String` also uses [`From`]; converting to
/// `String` uses [`TryFrom`] to check UTF-8 in each field.
/// `Box<[u8]>` fields also convert to `String` with [`TryFrom`], reusing their allocations.
/// Borrowed `&str` fields use the target's native [`From`] implementation:
/// `String`, `Vec<u8>` and [`BytesMut`] copy the text, while [`Bytes`] borrows
/// `&'static str` fields without copying.
/// Borrowed `&[u8]` fields likewise copy into `Vec<u8>`, [`BytesMut`] or `Box<[u8]>`;
/// `&'static [u8]` fields convert to [`Bytes`] without copying.
/// Converting `&[u8]` fields to `String` uses [`TryFrom`]: each field is copied
/// into owned storage and checked for UTF-8.
/// These conversions also support [`TryInto`]; infallible conversions use
/// [`std::convert::Infallible`] as their error type.
/// To borrow byte fields as `&str`, use `(&fields).try_into()`. This checks
/// UTF-8 without copying the field bytes, and the result borrows from `fields`.
/// This supports `&[u8]`, `Vec<u8>`, `Box<[u8]>`, [`Bytes`] and [`BytesMut`].
/// For other borrowed fields, first use `FastqField::<&[u8]>::from(&fields)`.
/// For `String` fields, `(&fields).into()` borrows `&str` fields without copying
/// or checking UTF-8 again. Its [`TryInto`] conversion has error type
/// [`std::convert::Infallible`].
/// Converting `FastqField<&[u8]>` to `FastqField<&str>` consumes the wrapper and
/// keeps the lifetime of its input slices, without copying their bytes.
///
/// ```
/// use bytes::Bytes;
/// use biorill_fastx::{FastqField, FastqParts};
///
/// let fields = FastqField::try_from(FastqParts::new(
///     String::from("read1"),
///     None,
///     String::from("AC"),
///     String::from(""),
///     String::from("12"),
/// ))?;
/// let bytes = FastqField::<Bytes>::from(fields);
/// let borrowed: FastqField<&str> = (&bytes).try_into()?;
/// assert_eq!(*borrowed.id(), "read1");
/// assert_eq!(borrowed.to_string(), "@read1\nAC\n+\n12\n");
/// let text: FastqField<String> = bytes.try_into()?;
/// assert_eq!(text.id(), "read1");
/// # Ok::<(), biorill_fastx::FastqFieldError>(())
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
        // Constructors and conversions only produce the supported storage types
        // or references to them. Their fixed Deref implementations preserve the
        // validated bytes. Deref alone does not guarantee this for custom types;
        // admitting custom storage would require revisiting this conversion.
        FastqField {
            id: self.id.deref(),
            description: self.description.as_deref(),
            sequence: self.sequence.deref(),
            separator: self.separator.deref(),
            quality: self.quality.deref(),
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
    /// Use `id` as the identifier, without its leading `@` or description.
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing the fields if `id` contains
    /// LF, a space, or a tab.
    pub fn use_id(&mut self, id: T) -> Result<(), FastqFieldError> {
        validate_line("identifier", id.as_bytes())?;
        validate_id_delimiters(id.as_bytes())?;
        self.id = id;
        Ok(())
    }

    /// Use the supplied description, or remove it with `None`.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::LineFeed`] without changing the fields if the
    /// description contains LF.
    pub fn use_description(&mut self, description: Option<T>) -> Result<(), FastqFieldError> {
        if let Some(description) = &description {
            validate_line("description", description.as_bytes())?;
        }
        self.description = description;
        Ok(())
    }

    /// Use the content after the separator line's leading `+`.
    /// An empty value produces a bare `+` line; supplied content is never trimmed.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::LineFeed`] without changing the fields if the
    /// separator content contains LF.
    pub fn use_separator(&mut self, separator: T) -> Result<(), FastqFieldError> {
        validate_line("separator", separator.as_bytes())?;
        self.separator = separator;
        Ok(())
    }

    /// Use a new sequence with the same byte length as the current quality scores.
    /// To change both lengths, use [`Self::use_sequence_and_quality`].
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing the fields if the sequence
    /// contains LF or its byte length differs from the current quality length.
    pub fn use_sequence(&mut self, sequence: T) -> Result<(), FastqFieldError> {
        validate_line("sequence", sequence.as_bytes())?;
        validate_lengths(sequence.as_bytes(), self.quality.as_bytes())?;
        self.sequence = sequence;
        Ok(())
    }

    /// Use new quality scores with the same byte length as the current sequence.
    /// To change both lengths, use [`Self::use_sequence_and_quality`].
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing the fields if the quality
    /// contains LF or its byte length differs from the current sequence length.
    pub fn use_quality(&mut self, quality: T) -> Result<(), FastqFieldError> {
        validate_line("quality", quality.as_bytes())?;
        validate_lengths(self.sequence.as_bytes(), quality.as_bytes())?;
        self.quality = quality;
        Ok(())
    }

    /// Use a new sequence and quality together. Empty values are accepted.
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] without changing either field if a replacement
    /// contains LF or the replacement lengths differ.
    pub fn use_sequence_and_quality(
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

impl<'a, T: StringBytes> From<&'a FastqField<T>> for FastqField<&'a [u8]> {
    /// Borrow the validated fields without copying or revalidating them.
    /// Updating the returned fields does not change the source fields.
    fn from(field: &'a FastqField<T>) -> Self {
        Self {
            id: field.id.as_bytes(),
            description: field.description.as_ref().map(StringBytes::as_bytes),
            sequence: field.sequence.as_bytes(),
            separator: field.separator.as_bytes(),
            quality: field.quality.as_bytes(),
        }
    }
}

impl<'a> From<&'a FastqField<String>> for FastqField<&'a str> {
    /// Borrow each string as text without copying or revalidating it.
    /// The result borrows from the source fields; updating it does not change them.
    fn from(fields: &'a FastqField<String>) -> Self {
        fields.as_deref()
    }
}

// A blanket From<FastqField<T>> for FastqField<U> where T: Into<U> would overlap
// the standard library's From<T> for T when the storage types are equal.
// It would also allow conversions that change the validated field bytes.
// Keep these implementations concrete so each conversion preserves those bytes.
macro_rules! impl_field_from {
    ($( $(#[$doc:meta])* $source:ty => $target:ty, $convert:expr; )+) => {
        $(
            impl From<FastqField<$source>> for FastqField<$target> {
                $(#[$doc])*
                fn from(fields: FastqField<$source>) -> Self {
                    Self {
                        id: ($convert)(fields.id),
                        description: fields.description.map($convert),
                        sequence: ($convert)(fields.sequence),
                        separator: ($convert)(fields.separator),
                        quality: ($convert)(fields.quality),
                    }
                }
            }
        )+
    };
}

// Borrowed inputs follow the target's native From implementations. Bytes only
// accepts static slices/text this way; extending it to shorter borrows would
// require an additional copying conversion.
impl_field_from! {
    /// Copy the borrowed byte fields into vectors using `Vec::from`.
    &[u8] => Vec<u8>, Vec::from;
    /// Copy the borrowed byte fields into boxed slices using `Box::from`.
    &[u8] => Box<[u8]>, Box::from;
    /// Borrow the static byte fields as byte buffers without copying their bytes.
    &'static [u8] => Bytes, Bytes::from;
    /// Copy the borrowed byte fields into mutable buffers using `BytesMut::from`.
    &[u8] => BytesMut, BytesMut::from;

    /// Consume the fields and convert their vectors to boxed slices.
    /// Removing spare capacity may reallocate their storage.
    Vec<u8> => Box<[u8]>, Vec::into_boxed_slice;
    /// Consume the fields and transfer their vector storage without copying their bytes.
    Vec<u8> => Bytes, Bytes::from;
    /// Consume the fields and transfer their vector storage without copying their bytes.
    Vec<u8> => BytesMut, |data: Vec<u8>| BytesMut::from(Bytes::from(data));

    /// Consume the fields and transfer their boxed storage without copying their bytes.
    Box<[u8]> => Vec<u8>, Vec::from;
    /// Consume the fields and transfer their boxed storage without copying their bytes.
    Box<[u8]> => Bytes, Bytes::from;
    /// Consume the fields and transfer their boxed storage without copying their bytes.
    Box<[u8]> => BytesMut, |data: Box<[u8]>| BytesMut::from(Bytes::from(data));

    /// Consume the fields and convert their bytes to vectors.
    /// Uniquely owned reusable storage is retained; shared or static bytes are copied.
    Bytes => Vec<u8>, Vec::from;
    /// Consume the fields and convert their bytes to boxed slices.
    /// Shared or static bytes are copied; removing spare capacity may reallocate storage.
    Bytes => Box<[u8]>, |data: Bytes| Vec::from(data).into_boxed_slice();
    /// Consume the fields and convert their bytes to mutable storage.
    /// Uniquely owned reusable storage is retained; shared or static bytes are copied.
    Bytes => BytesMut, BytesMut::from;

    /// Consume the fields and convert their buffers to vectors.
    /// Uniquely owned storage is reused; shared storage is copied.
    BytesMut => Vec<u8>, Vec::from;
    /// Consume the fields and convert their buffers to boxed slices.
    /// Shared storage is copied; removing spare capacity may reallocate storage.
    BytesMut => Box<[u8]>, |data: BytesMut| Vec::from(data).into_boxed_slice();
    /// Consume the fields and freeze their buffers without copying their bytes.
    BytesMut => Bytes, BytesMut::freeze;

    /// Copy the borrowed text fields into byte vectors using `Vec::from`.
    &str => Vec<u8>, Vec::from;
    /// Borrow the static text fields as byte buffers without copying their bytes.
    &'static str => Bytes, Bytes::from;
    /// Copy the borrowed text fields into mutable buffers using `BytesMut::from`.
    &str => BytesMut, BytesMut::from;
    /// Copy the borrowed text fields into owned strings using `String::from`.
    &str => String, String::from;

    /// Consume the fields and transfer their string storage without copying their bytes.
    String => Vec<u8>, String::into_bytes;
    /// Consume the fields and convert their strings to boxed byte slices.
    /// Removing spare capacity may reallocate their storage.
    String => Box<[u8]>, |data: String| data.into_bytes().into_boxed_slice();
    /// Consume the fields and transfer their string storage without copying their bytes.
    String => Bytes, Bytes::from;
    /// Consume the fields and transfer their string storage without copying their bytes.
    String => BytesMut, |data: String| BytesMut::from(Bytes::from(data));
}

// Keep these implementations concrete: the text From conversions also supply
// TryFrom through the standard library, which overlaps a generic StringBytes impl.
macro_rules! impl_field_try_as_str {
    ($($source:ty),+ $(,)?) => {
        $(
            impl<'a> TryFrom<&'a FastqField<$source>> for FastqField<&'a str> {
                type Error = FastqFieldError;

                /// Check UTF-8 and borrow each field as text without copying it.
                /// The result borrows from the source fields.
                ///
                /// # Errors
                /// Returns [`FastqFieldError::InvalidUtf8Slice`] for the first field
                /// containing invalid UTF-8.
                fn try_from(fields: &'a FastqField<$source>) -> Result<Self, Self::Error> {
                    fields.as_deref().try_into()
                }
            }
        )+
    };
}

impl_field_try_as_str!(&[u8], Vec<u8>, Box<[u8]>, Bytes, BytesMut);

// This by-value conversion consumes only slice references; their backing bytes
// outlive the wrapper. A generic by-value conversion could instead drop owned
// storage before the returned text is used.
impl<'a> TryFrom<FastqField<&'a [u8]>> for FastqField<&'a str> {
    type Error = FastqFieldError;

    /// Consume the borrowed fields and check each slice for UTF-8.
    /// The text keeps the input slices' lifetime without copying their bytes.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::InvalidUtf8Slice`] for the first field
    /// containing invalid UTF-8.
    fn try_from(fields: FastqField<&'a [u8]>) -> Result<Self, Self::Error> {
        let from_utf8 = |field, bytes| {
            StringBytes::try_borrow_str(bytes)
                .map_err(|source| FastqFieldError::InvalidUtf8Slice { field, source })
        };
        Ok(Self {
            id: from_utf8("identifier", fields.id)?,
            description: fields
                .description
                .map(|bytes| from_utf8("description", bytes))
                .transpose()?,
            sequence: from_utf8("sequence", fields.sequence)?,
            separator: from_utf8("separator", fields.separator)?,
            quality: from_utf8("quality", fields.quality)?,
        })
    }
}

// StringBytes preserves the bytes during text conversion. Keep concrete
// sources because a generic TryFrom would overlap existing From conversions.
macro_rules! impl_field_try_into_string {
    ($($source:ty),+ $(,)?) => {
        $(
            impl TryFrom<FastqField<$source>> for FastqField<String> {
                type Error = FastqFieldError;

                /// Consume the fields and convert each component to UTF-8 text.
                /// Borrowed slices are copied. An absent description stays absent.
                ///
                /// # Errors
                /// Returns [`FastqFieldError::InvalidUtf8`] for the first field
                /// containing invalid UTF-8, retaining its name and rejected bytes.
                fn try_from(fields: FastqField<$source>) -> Result<Self, Self::Error> {
                    let from_utf8 = |field, data: $source| {
                        data.try_into_string()
                            .map_err(|source| FastqFieldError::InvalidUtf8 { field, source })
                    };
                    Ok(Self {
                        id: from_utf8("identifier", fields.id)?,
                        description: fields.description
                            .map(|data| from_utf8("description", data))
                            .transpose()?,
                        sequence: from_utf8("sequence", fields.sequence)?,
                        separator: from_utf8("separator", fields.separator)?,
                        quality: from_utf8("quality", fields.quality)?,
                    })
                }
            }
        )+
    };
}

impl_field_try_into_string!(&[u8], Vec<u8>, Box<[u8]>, Bytes, BytesMut);

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
