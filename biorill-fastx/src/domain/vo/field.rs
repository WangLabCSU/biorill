use std::fmt;
use std::io::{self, Write};
use std::ops::Deref;
use std::str::Utf8Error;
use std::string::FromUtf8Error;

#[cfg(feature = "bytes")]
use bytes::BufMut;
use memchr::{memchr, memchr2};
use stringbytes::StringBytes;
use thiserror::Error;

use super::parts::FastqParts;
use super::record::FastqRecord;

/// A FASTQ record with individually editable fields.
///
/// Create a record from [`FastqParts`] with [`Self::try_from`], then use the
/// `set_*` methods to replace fields. Construction and updates check the FASTQ
/// field constraints, including equal sequence and quality byte lengths.
/// Use [`Self::check_sequence`] and [`Self::check_quality`] for optional
/// nucleotide and quality-character checks.
///
/// Fields can use any storage supported by [`StringBytes`].
///
/// # Examples
///
/// ```
/// use biorill_fastx::{FastqField, FastqParts};
///
/// let mut fields = FastqField::try_from(FastqParts::new("read1", None, "AC", "", "12"))?;
/// fields.set_description(Some("UMI:AC"))?;
/// fields.set_sequence_and_quality("C", "2")?;
/// assert_eq!(fields.to_string(), "@read1 UMI:AC\nC\n+\n2\n");
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
    /// Returns the identifier, excluding the leading `@` and description.
    pub fn id(&self) -> &T {
        &self.id
    }

    /// Returns the description, excluding the whitespace that separates it from the identifier.
    pub fn description(&self) -> Option<&T> {
        self.description.as_ref()
    }

    /// Returns the sequence.
    pub fn sequence(&self) -> &T {
        &self.sequence
    }

    /// Returns the content after the separator line's `+`. A bare `+` line has empty content.
    pub fn separator(&self) -> &T {
        &self.separator
    }

    /// Returns the quality scores, one byte per sequence byte.
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

    /// Borrows each field, returning `FastqField<&T>`.
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
    /// assert_eq!(*borrowed.id(), "read1");
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

    /// Returns the fields as editable [`FastqParts`], consuming `self`.
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
    /// Checks the sequence against the IUPAC nucleotide alphabet `ACGTURYSWKMBDHVN`.
    /// Accepts uppercase and lowercase letters. Empty sequences are accepted.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::InvalidSequence`] with the first invalid byte and
    /// its zero-based position in the sequence.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new("r1", None, "AcUn", "", "!I~I"))?;
    /// fields.check_sequence()?;
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn check_sequence(&self) -> Result<(), FastqFieldError> {
        // Bits 0 through 25 correspond to uppercase A through Z.
        const IUPAC_NUCLEOBASE_BITSET: u32 = {
            let alphabet = b"ACGTURYSWKMBDHVN";
            let mut bitset = 0u32;
            let mut i = 0;
            while i < alphabet.len() {
                // Normalize ASCII lowercase letters to uppercase, then map 'A'..='Z'
                // to indices 0..=25. Wrapping subtraction maps bytes below 'A' to
                // values outside this range, allowing the bounds check to reject them.
                bitset |= 1u32 << (alphabet[i] - b'A');
                i += 1;
            }
            bitset
        };

        for (index, &byte) in self.sequence.as_bytes().iter().enumerate() {
            let pos = byte.to_ascii_uppercase().wrapping_sub(b'A');
            if pos >= 26 || IUPAC_NUCLEOBASE_BITSET & (1u32 << pos) == 0 {
                return Err(FastqFieldError::InvalidSequence { index, byte });
            }
        }
        Ok(())
    }

    /// Checks that quality bytes are in the Phred+33 range `!` through `~` (33–126).
    /// Scans the quality scores without allocating. Empty quality scores are accepted.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::InvalidQuality`] with the first invalid byte and
    /// its zero-based position in the quality scores.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqFieldError, FastqParts};
    ///
    /// let mut fields = FastqField::try_from(FastqParts::new("r1", None, "ACU", "", "!I~"))?;
    /// fields.check_quality()?;
    ///
    /// fields.set_quality("! ~")?;
    /// assert!(matches!(
    ///     fields.check_quality(),
    ///     Err(FastqFieldError::InvalidQuality { index: 1, byte: b' ' })
    /// ));
    /// # Ok::<(), FastqFieldError>(())
    /// ```
    pub fn check_quality(&self) -> Result<(), FastqFieldError> {
        // SWAR checks several bytes with ordinary integer operations: one usize
        // holds 8 bytes on a 64-bit target, or 4 on a 32-bit target.
        const WORD_BYTES: usize = (usize::BITS / 8) as usize;
        // Repeat each constant in every byte of the word. b'\x01' is the value
        // 1; b'\x80' is 0b1000_0000, selecting each byte's highest bit.
        const LOWER_BOUND: usize = usize::from_ne_bytes([b'!'; WORD_BYTES]);
        const INCREMENT: usize = usize::from_ne_bytes([b'\x01'; WORD_BYTES]);
        const HIGH_BITS: usize = usize::from_ne_bytes([b'\x80'; WORD_BYTES]);

        let quality = self.quality.as_bytes();
        let (chunks, tail) = quality.as_chunks::<WORD_BYTES>();
        let invalid_chunk = chunks.iter().position(|chunk| {
            // Interpret this byte array as one native-endian integer. This packs
            // the bytes for arithmetic; it does not parse their text as a number.
            let word = usize::from_ne_bytes(*chunk);
            // For a byte b with no incoming borrow/carry, the high bits are:
            //
            // b range                 | b | b - b'!' | b + b'\x01'
            // ------------------------+---+----------+-------------
            // 0..=32 (below b'!')      | 0 |    1     |     0
            // b'!'..=b'~' (33..=126)   | 0 |    0     |     0
            // b'\x7f' (127)           | 0 |    0     |     1
            // b'\x80'..=b'\xff'       | 1 |    *     |     *
            //
            // Subtraction wraps bytes below b'!' into 223..=255; adding 1
            // turns 127 into 128. Bytes >=128 already have their high bit set.
            // OR combines these three checks; * means that value is irrelevant.
            //
            // Arithmetic acts on the whole word, so borrows/carries can cross
            // byte boundaries. All-valid bytes produce neither. Otherwise, the
            // least significant invalid byte has only valid bytes below it,
            // so it receives no borrow/carry and is detected by the table.
            // Higher bytes can also be flagged, including valid ones. The mask
            // therefore detects a bad block, not the exact set of bad bytes.
            // This any-error test works with either native byte order.
            let invalid = word | word.wrapping_sub(LOWER_BOUND) | word.wrapping_add(INCREMENT);
            invalid & HIGH_BITS != 0
        });

        // position() returns a block index, not the invalid byte's offset.
        // For example, b"III IIII" fails in block 0, but its bad byte is at 3.
        // Scan that block in input order to find the first error, stopping
        // within WORD_BYTES bytes. Earlier blocks are not scanned again.
        // If all full blocks passed, only the remaining short tail is checked.
        let offset = invalid_chunk.map_or(quality.len() - tail.len(), |index| index * WORD_BYTES);
        for (index, &byte) in quality[offset..].iter().enumerate() {
            if !byte.is_ascii_graphic() {
                return Err(FastqFieldError::InvalidQuality {
                    index: offset + index,
                    byte,
                });
            }
        }
        Ok(())
    }

    /// Borrows each field through [`Deref`].
    /// For example, `FastqField<Vec<u8>>` becomes `FastqField<&[u8]>`.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     b"read1".to_vec(), None, b"AC".to_vec(), Vec::new(), b"12".to_vec(),
    /// ))?;
    /// let slices: FastqField<&[u8]> = fields.as_deref();
    /// assert_eq!(*slices.sequence(), b"AC");
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

    /// Borrows all fields as byte slices.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     "read1", None, "AC", "", "12",
    /// ))?;
    /// let view = fields.byte_view();
    /// assert_eq!(*view.sequence(), b"AC");
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn byte_view(&self) -> FastqField<&[u8]> {
        // SAFETY: StringBytes::as_bytes exposes each component's validated bytes.
        unsafe { self.as_ref().map_fields(T::as_bytes) }
    }

    /// Borrows all fields as UTF-8 text.
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

    /// Replaces the identifier. Supply only the name, without `@` or a description.
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

    /// Replaces the description, or removes it with `None`.
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

    /// Replaces the content after the separator line's `+`.
    /// Use an empty value for a bare `+` line. The supplied content is kept as-is.
    ///
    /// # Errors
    /// Returns [`FastqFieldError::LineFeed`] without changing the fields if the
    /// separator content contains LF.
    pub fn set_separator(&mut self, separator: T) -> Result<(), FastqFieldError> {
        validate_line("separator", separator.as_bytes())?;
        self.separator = separator;
        Ok(())
    }

    /// Replaces the sequence, keeping its current byte length.
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

    /// Replaces the quality scores, keeping their current byte length.
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

    /// Replaces both the sequence and quality scores, allowing their length to change.
    /// The new values must have equal byte lengths; both may be empty.
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

    /// Writes the record as four LF-terminated FASTQ lines.
    /// A present description follows the identifier, separated by one space.
    /// Field contents are written as-is, including non-UTF-8 bytes.
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

    /// Appends the record to `buf` in the same FASTQ format as [`Self::write`].
    /// Copies the field contents and grows the vector if needed.
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
        buf.push(b'@');
        buf.extend_from_slice(self.id.as_bytes());
        if let Some(description) = &self.description {
            buf.push(b' ');
            buf.extend_from_slice(description.as_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(self.sequence.as_bytes());
        buf.extend_from_slice(b"\n+");
        buf.extend_from_slice(self.separator.as_bytes());
        buf.push(b'\n');
        buf.extend_from_slice(self.quality.as_bytes());
        buf.push(b'\n');
    }

    /// Appends the record to `buf` in the same FASTQ format as [`Self::write`].
    /// Copies the field contents into the buffer.
    ///
    /// Available with the `bytes` feature.
    ///
    /// # Panics
    /// Panics if `buf` cannot accept the fields. A prefix may already have been appended.
    #[cfg(feature = "bytes")]
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
    /// Converts the fields to owned byte vectors.
    /// Consumes `fields`, reusing storage where possible and copying otherwise.
    /// See [`StringBytes::into_bytes`] for the costs of each storage type.
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts};
    ///
    /// let fields = FastqField::try_from(FastqParts::new(
    ///     String::from("r1"), None, String::from("AC"),
    ///     String::new(), String::from("12"),
    /// ))?;
    /// let bytes = FastqField::<Vec<u8>>::from_fields(fields);
    /// assert_eq!(bytes.sequence(), b"AC");
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    pub fn from_fields<T: StringBytes>(fields: FastqField<T>) -> Self {
        // SAFETY: StringBytes::into_bytes preserves each component's bytes.
        unsafe { fields.map_fields(StringBytes::into_bytes) }
    }
}

impl FastqField<String> {
    /// Converts the fields to owned strings, checking UTF-8 where needed.
    /// Consumes `fields`, reusing storage where possible and copying otherwise.
    /// See [`StringBytes::try_into_string`] for the costs of each storage type.
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
    /// Borrows the record's fields as byte slices.
    /// Replacing fields in the result does not change the source record.
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
    /// Formats the fields as four LF-terminated FASTQ lines.
    /// Uses each field's [`fmt::Display`] implementation, with a space before any description.
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
    /// The sequence contains a byte outside the accepted nucleotide alphabet.
    #[error("invalid FASTQ sequence byte {byte:#04x} at position {index}")]
    InvalidSequence {
        /// Zero-based byte position in the sequence.
        index: usize,
        /// The rejected byte.
        byte: u8,
    },
    /// A quality byte is outside the Phred+33 range `!` through `~` (33–126).
    #[error("invalid FASTQ quality byte {byte:#04x} at position {index}")]
    InvalidQuality {
        /// Zero-based byte position in the quality scores.
        index: usize,
        /// The rejected byte.
        byte: u8,
    },
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

    /// Creates a record by checking the supplied fields.
    /// Scans their contents and takes ownership without copying them.
    ///
    /// Supply `separator` without the line's leading `+`; an empty value gives
    /// a bare `+` line. Values are kept as-is, including any additional `+`.
    /// Empty identifiers and empty sequence/quality pairs are accepted.
    /// To also check nucleotide and quality characters, call [`Self::check_sequence`]
    /// and [`Self::check_quality`] on the result.
    ///
    /// # Errors
    /// Returns [`FastqFieldError`] if a field contains LF, the identifier contains
    /// a space or tab, or the sequence and quality have different byte lengths.
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
    if memchr(b'\n', bytes).is_some() {
        return Err(FastqFieldError::LineFeed { field });
    }
    Ok(())
}

fn validate_id_delimiters(id: &[u8]) -> Result<(), FastqFieldError> {
    if memchr2(b' ', b'\t', id).is_some() {
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

    #[cfg(feature = "bytes")]
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
    fn optional_checks_accept_only_the_documented_characters() {
        for byte in 0..=u8::MAX {
            if byte == b'\n' {
                continue; // Construction already rejects line feeds.
            }
            let input = [byte];
            let fields = FastqField::try_from(FastqParts::new(
                b"r".as_slice(),
                None,
                input.as_slice(),
                b"",
                input.as_slice(),
            ))
            .unwrap();

            let sequence = fields.check_sequence();
            if b"ACGTURYSWKMBDHVNacgturyswkmbdhvn".contains(&byte) {
                assert!(sequence.is_ok(), "sequence byte {byte:#04x}");
            } else {
                assert!(matches!(
                    sequence,
                    Err(FastqFieldError::InvalidSequence { index: 0, byte: rejected })
                        if rejected == byte
                ));
            }

            let quality = fields.check_quality();
            if byte.is_ascii_graphic() {
                assert!(quality.is_ok(), "quality byte {byte:#04x}");
            } else {
                assert!(matches!(
                    quality,
                    Err(FastqFieldError::InvalidQuality { index: 0, byte: rejected })
                        if rejected == byte
                ));
            }
        }
        let empty = FastqField::try_from(FastqParts::new("r", None, "", "", "")).unwrap();
        assert!(empty.check_sequence().is_ok());
        assert!(empty.check_quality().is_ok());
    }

    #[test]
    fn optional_checks_report_the_first_invalid_byte_in_the_requested_field() {
        let mut fields = FastqField::try_from(FastqParts::new(
            b"\xff".as_slice(),
            None,
            b"ACUN",
            b"\xff",
            b"!I~I",
        ))
        .unwrap();
        fields.set_sequence(b"AU?X").unwrap();
        assert!(fields.check_quality().is_ok());
        let error = fields.check_sequence().unwrap_err();
        assert!(matches!(
            error,
            FastqFieldError::InvalidSequence {
                index: 2,
                byte: b'?'
            }
        ));
        assert_eq!(
            error.to_string(),
            "invalid FASTQ sequence byte 0x3f at position 2"
        );

        fields.set_sequence(b"ACUN").unwrap();
        fields.set_quality(b"! \x7f~").unwrap();
        assert!(fields.check_sequence().is_ok());
        let error = fields.check_quality().unwrap_err();
        assert!(matches!(
            error,
            FastqFieldError::InvalidQuality {
                index: 1,
                byte: b' '
            }
        ));
        assert_eq!(
            error.to_string(),
            "invalid FASTQ quality byte 0x20 at position 1"
        );
    }

    #[test]
    fn sequence_checks_locate_errors_across_block_boundaries_and_tails() {
        let alphabet = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn";
        for length in [0, 1, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 500] {
            let mut sequence: Vec<u8> = alphabet.iter().copied().cycle().take(length).collect();
            let quality = vec![b'!'; length];
            let check = |sequence: &[u8]| {
                FastqField::try_from(FastqParts::new(
                    b"r".as_slice(),
                    None,
                    sequence,
                    b"",
                    quality.as_slice(),
                ))
                .unwrap()
                .check_sequence()
            };
            assert!(check(&sequence).is_ok());
            for position in 0..length {
                sequence[length - 1] = b'x';
                sequence[position] = b'?';
                assert!(matches!(
                    check(&sequence),
                    Err(FastqFieldError::InvalidSequence { index, byte: b'?' })
                        if index == position
                ));
                sequence[position] = alphabet[position % alphabet.len()];
                sequence[length - 1] = alphabet[(length - 1) % alphabet.len()];
            }
        }
    }

    #[test]
    fn quality_checks_locate_errors_across_block_boundaries_and_tails() {
        for length in [0, 1, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65] {
            let sequence = vec![b'A'; length];
            let mut quality: Vec<u8> = (b'!'..=b'~').cycle().take(length).collect();
            let check = |quality: &[u8]| {
                FastqField::try_from(FastqParts::new(
                    b"r".as_slice(),
                    None,
                    sequence.as_slice(),
                    b"",
                    quality,
                ))
                .unwrap()
                .check_quality()
            };
            assert!(check(&quality).is_ok());
            for position in 0..length {
                let original = quality[position];
                for byte in 0..=u8::MAX {
                    if byte == b'\n' {
                        continue; // Construction already rejects line feeds.
                    }
                    quality[position] = byte;
                    if (33..=126).contains(&byte) {
                        assert!(check(&quality).is_ok());
                    } else {
                        assert!(matches!(
                            check(&quality),
                            Err(FastqFieldError::InvalidQuality { index, byte: rejected })
                                if index == position && rejected == byte
                        ));
                    }
                }
                quality[position] = original;
            }
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
        #[cfg(feature = "bytes")]
        check(Bytes::from);
        #[cfg(feature = "bytes")]
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
        #[cfg(feature = "bytes")]
        check(Bytes::from_static);
        #[cfg(feature = "bytes")]
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

    #[cfg(feature = "bytes")]
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
        #[cfg(feature = "bytes")]
        check(|text| Bytes::copy_from_slice(text.as_bytes()));
        #[cfg(feature = "bytes")]
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
        #[cfg(feature = "bytes")]
        check(Bytes::from);
        #[cfg(feature = "bytes")]
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
        #[cfg(feature = "bytes")]
        check(Bytes::from);
        #[cfg(feature = "bytes")]
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
            record.encode_into_vec(&mut encoded);
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
        let original = FastqRecord::<Vec<u8>>::try_from(valid_parts()).unwrap();
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
        let updated = FastqRecord::<Vec<u8>>::from(fields);
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
            FastqRecord::<&[u8]>::from(&FastqRecord::<Vec<u8>>::from(fields))
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
        let updated = FastqRecord::<Vec<u8>>::from(FastqField::try_from(parts).unwrap());
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
