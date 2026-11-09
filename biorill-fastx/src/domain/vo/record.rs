use std::io::{self, Write};
use std::ops::Range;
use std::str::Utf8Error;
use std::string::FromUtf8Error;

use bytes::{BufMut, Bytes, BytesMut};
use stringbytes::StringBytes;
use thiserror::Error;

use super::field::{FastqField, FastqFieldError};
use super::parts::FastqParts;

/// A FASTQ record with an identifier, optional description, sequence, separator
/// and quality scores. Sequence and quality have equal byte lengths.
///
/// To display a record as FASTQ text, borrow its fields with [`Self::fields`]
/// and check UTF-8 with [`FastqField::try_str_view`] to obtain `FastqField<&str>`.
/// Field ranges are byte ranges: even text storage can select only part of a
/// UTF-8 character, so its fields still need this check before text output.
///
/// Construct from a buffer and field ranges with [`Self::try_from`]. The buffer
/// must implement [`StringBytes`], which provides a consistent byte view for
/// validated ranges. Byte and text storage, references, shared pointers, and
/// nested wrappers are supported. Construction checks the selected fields and
/// retains the buffer without copying it.
///
/// ```
/// use std::rc::Rc;
/// use biorill_fastx::{FastqParts, FastqRecord};
///
/// let record = FastqRecord::try_from((
///     Rc::new(String::from("read1AC12")),
///     FastqParts::new(0..5, None, 5..7, 7..7, 7..9),
/// ))?;
/// assert_eq!(record.sequence(), b"AC");
/// assert_eq!(record.quality(), b"12");
/// # Ok::<(), biorill_fastx::FastqRecordError>(())
/// ```
///
/// Custom types implementing only `AsRef<[u8]>` must lend a byte slice for
/// validation instead of being retained as the record's buffer.
///
/// ```compile_fail,E0277
/// use biorill_fastx::{FastqParts, FastqRecord};
///
/// struct Custom(Vec<u8>);
/// impl AsRef<[u8]> for Custom {
///     fn as_ref(&self) -> &[u8] { &self.0 }
/// }
///
/// let record = FastqRecord::try_from((
///     Custom(b"read1AC12".to_vec()),
///     FastqParts::new(0..5, None, 5..7, 7..7, 7..9),
/// ));
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqRecord<T> {
    buf: T,
    ranges: FastqParts<Range<usize>>,
}

impl<T> FastqRecord<T> {
    /// Borrow the backing buffer, converting `&FastqRecord<T>` to `FastqRecord<&T>`.
    /// Retains the field ranges without copying or revalidating the bytes.
    ///
    /// ```
    /// use biorill_fastx::{FastqParts, FastqRecord};
    ///
    /// let record = FastqRecord::try_from((
    ///     String::from("read1AC12"),
    ///     FastqParts::new(0..5, None, 5..7, 7..7, 7..9),
    /// ))?;
    /// let borrowed: FastqRecord<&String> = record.as_ref();
    /// assert_eq!(borrowed.sequence(), b"AC");
    /// assert_eq!(borrowed.into_parts().0, "read1AC12");
    /// # Ok::<(), biorill_fastx::FastqRecordError>(())
    /// ```
    pub fn as_ref(&self) -> FastqRecord<&T> {
        FastqRecord {
            buf: &self.buf,
            ranges: self.ranges.clone(),
        }
    }

    /// Consume the record and return its backing buffer and field ranges.
    pub fn into_parts(self) -> (T, FastqParts<Range<usize>>) {
        (self.buf, self.ranges)
    }
}

impl<T: StringBytes> FastqRecord<T> {
    /// The identifier without the leading `@` or description.
    pub fn id(&self) -> &[u8] {
        &self.buf.as_bytes()[self.ranges.id.clone()]
    }

    /// The optional description without the identifier or separating whitespace.
    pub fn description(&self) -> Option<&[u8]> {
        self.ranges
            .description
            .as_ref()
            .map(|range| &self.buf.as_bytes()[range.clone()])
    }

    /// The sequence in input order.
    pub fn sequence(&self) -> &[u8] {
        &self.buf.as_bytes()[self.ranges.sequence.clone()]
    }

    /// The content after the separator line's leading `+`, possibly empty.
    pub fn separator(&self) -> &[u8] {
        &self.buf.as_bytes()[self.ranges.separator.clone()]
    }

    /// The quality scores in the same order as the sequence.
    pub fn quality(&self) -> &[u8] {
        &self.buf.as_bytes()[self.ranges.quality.clone()]
    }

    /// Borrow this record's validated fields without copying or revalidating them.
    /// Updating the returned fields does not change this record.
    ///
    /// To display the fields as FASTQ text, explicitly check UTF-8 first:
    ///
    /// ```
    /// use biorill_fastx::{FastqField, FastqParts, FastqRecord};
    ///
    /// let record = FastqRecord::try_from((
    ///     b"read1AC12".as_slice(),
    ///     FastqParts::new(0..5, None, 5..7, 7..7, 7..9),
    /// ))?;
    /// let fields = record.fields();
    /// let text: FastqField<&str> = fields.try_str_view()?;
    /// assert_eq!(text.to_string(), "@read1\nAC\n+\n12\n");
    /// # Ok::<(), biorill_fastx::FastqRecordError>(())
    /// ```
    pub fn fields(&self) -> FastqField<&[u8]> {
        FastqField::from(self)
    }

    /// Borrow the components for construction of another record.
    /// Changes to the returned parts do not modify this record.
    pub fn parts(&self) -> FastqParts<&[u8]> {
        FastqParts {
            id: self.id(),
            description: self.description(),
            sequence: self.sequence(),
            separator: self.separator(),
            quality: self.quality(),
        }
    }

    /// Write four FASTQ lines with LF endings and a space before any description.
    /// Field bytes are written without checking their encoding.
    ///
    /// # Errors
    /// Returns the writer's error if the record cannot be written in full.
    /// A prefix of the record may already have been written.
    pub fn write<W: Write + ?Sized>(&self, writer: &mut W) -> io::Result<()> {
        self.fields().write(writer)
    }

    /// Encode four FASTQ lines and append them to `buf`, using LF endings and a
    /// space before any description.
    /// Field bytes are copied into `buf` without checking their encoding.
    ///
    /// # Panics
    /// Panics if `buf` cannot accept the record. A prefix may already have been appended.
    pub fn encode_into_buf<B: BufMut + ?Sized>(&self, buf: &mut B) {
        self.fields().encode_into_buf(buf);
    }

    /// Encode this record and append it to `buf`, reserving space for the
    /// complete record first. Uses LF endings and a space before any description.
    /// Copies all field bytes; growing `buf` may allocate.
    pub fn encode_into_vec(&self, buf: &mut Vec<u8>) {
        self.fields().encode_into_vec(buf);
    }
}

impl<'a, T: StringBytes> From<&'a FastqRecord<T>> for FastqRecord<&'a [u8]> {
    /// Borrow the backing buffer and field ranges without copying the bytes or
    /// revalidating the record.
    fn from(record: &'a FastqRecord<T>) -> Self {
        Self {
            buf: record.buf.as_bytes(),
            ranges: record.ranges.clone(),
        }
    }
}

// A blanket From<FastqRecord<T>> for FastqRecord<U> where T: Into<U> would overlap
// the standard library's From<T> for T. Arbitrary Into implementations could
// also change the bytes and invalidate the stored ranges. These concrete
// conversions preserve the complete byte contents and order.
macro_rules! impl_record_from {
    ($( $(#[$doc:meta])* $source:ty => $target:ty, $convert:expr; )+) => {
        $(
            impl From<FastqRecord<$source>> for FastqRecord<$target> {
                $(#[$doc])*
                fn from(record: FastqRecord<$source>) -> Self {
                    Self {
                        buf: ($convert)(record.buf),
                        ranges: record.ranges,
                    }
                }
            }
        )+
    };
}

// Borrowed inputs use native From implementations. Bytes accepts only static
// slices/text without an explicit copying operation.
impl_record_from! {
    /// Copy the borrowed storage into a vector using `Vec::from`.
    &[u8] => Vec<u8>, Vec::from;
    /// Copy the borrowed storage into a boxed slice using `Box::from`.
    &[u8] => Box<[u8]>, Box::from;
    /// Borrow the static storage as a byte buffer without copying its bytes.
    &'static [u8] => Bytes, Bytes::from;
    /// Copy the borrowed storage into a mutable buffer using `BytesMut::from`.
    &[u8] => BytesMut, BytesMut::from;

    /// Consume the record and convert its vector to a boxed slice.
    /// Removing spare capacity may reallocate its storage.
    Vec<u8> => Box<[u8]>, Vec::into_boxed_slice;
    /// Consume the record and transfer its vector storage without copying its bytes.
    Vec<u8> => Bytes, Bytes::from;
    /// Consume the record and transfer its vector storage without copying its bytes.
    Vec<u8> => BytesMut, |data: Vec<u8>| BytesMut::from(Bytes::from(data));

    /// Consume the record and transfer its boxed storage without copying its bytes.
    Box<[u8]> => Vec<u8>, Vec::from;
    /// Consume the record and transfer its boxed storage without copying its bytes.
    Box<[u8]> => Bytes, Bytes::from;
    /// Consume the record and transfer its boxed storage without copying its bytes.
    Box<[u8]> => BytesMut, |data: Box<[u8]>| BytesMut::from(Bytes::from(data));

    /// Consume the record and convert its bytes to a vector.
    /// Reuses uniquely owned reusable storage, moving bytes to its beginning if needed.
    /// Shared or static bytes are copied.
    Bytes => Vec<u8>, Vec::from;
    /// Consume the record and convert its bytes to a boxed slice.
    /// Shared or static bytes are copied; removing spare capacity may reallocate storage.
    Bytes => Box<[u8]>, |data: Bytes| Vec::from(data).into_boxed_slice();
    /// Consume the record and convert its bytes to mutable storage.
    /// Uniquely owned reusable storage is retained; shared or static bytes are copied.
    Bytes => BytesMut, BytesMut::from;

    /// Consume the record and convert its bytes to a vector.
    /// Reuses uniquely owned storage, moving bytes to its beginning if needed.
    /// Shared storage is copied.
    BytesMut => Vec<u8>, Vec::from;
    /// Consume the record and convert its buffer to a boxed slice.
    /// Shared storage is copied; removing spare capacity may reallocate storage.
    BytesMut => Box<[u8]>, |data: BytesMut| Vec::from(data).into_boxed_slice();
    /// Consume the record and freeze its buffer without copying its bytes.
    BytesMut => Bytes, BytesMut::freeze;

    /// Copy the borrowed text into a byte vector using `Vec::from`.
    &str => Vec<u8>, Vec::from;
    /// Borrow the static text as a byte buffer without copying its bytes.
    &'static str => Bytes, Bytes::from;
    /// Copy the borrowed text into a mutable buffer using `BytesMut::from`.
    &str => BytesMut, BytesMut::from;
    /// Copy the borrowed text into an owned string using `String::from`.
    &str => String, String::from;

    /// Consume the record and transfer its string storage without copying its bytes.
    String => Vec<u8>, String::into_bytes;
    /// Consume the record and convert its string to a boxed byte slice.
    /// Removing spare capacity may reallocate its storage.
    String => Box<[u8]>, |data: String| data.into_bytes().into_boxed_slice();
    /// Consume the record and transfer its string storage without copying its bytes.
    String => Bytes, Bytes::from;
    /// Consume the record and transfer its string storage without copying its bytes.
    String => BytesMut, |data: String| BytesMut::from(Bytes::from(data));
}

// StringBytes text conversion preserves the buffer bytes and therefore its
// ranges. A concrete From for borrowed text records would overlap this TryFrom.
impl<'a, T: StringBytes> TryFrom<&'a FastqRecord<T>> for FastqRecord<&'a str> {
    type Error = FastqRecordError;

    /// Check UTF-8 and borrow the complete backing storage as text without copying it.
    /// Field ranges remain byte ranges and need not lie on character boundaries.
    ///
    /// ```
    /// use biorill_fastx::{FastqParts, FastqRecord};
    ///
    /// let record = FastqRecord::try_from((
    ///     b"read1AC12".to_vec(),
    ///     FastqParts::new(0..5, None, 5..7, 7..7, 7..9),
    /// ))?;
    /// let text: FastqRecord<&str> = (&record).try_into()?;
    /// assert_eq!(text.into_parts().0, "read1AC12");
    /// # Ok::<(), biorill_fastx::FastqRecordError>(())
    /// ```
    ///
    /// # Errors
    /// Returns [`FastqRecordError::InvalidUtf8Slice`] if any backing bytes are
    /// invalid UTF-8, including bytes outside the field ranges.
    fn try_from(record: &'a FastqRecord<T>) -> Result<Self, Self::Error> {
        Ok(Self {
            buf: record.buf.try_borrow_str()?,
            ranges: record.ranges.clone(),
        })
    }
}

// Consuming a slice view preserves its backing storage's lifetime. Consuming
// arbitrary owned storage could instead drop the bytes before the text is used.
impl<'a> TryFrom<FastqRecord<&'a [u8]>> for FastqRecord<&'a str> {
    type Error = FastqRecordError;

    /// Consume the borrowed record and check its complete backing slice for UTF-8.
    /// The text keeps the input slice's lifetime without copying its bytes.
    /// Field ranges remain byte ranges and need not lie on character boundaries.
    ///
    /// # Errors
    /// Returns [`FastqRecordError::InvalidUtf8Slice`] if any backing bytes are
    /// invalid UTF-8, including bytes outside the field ranges.
    fn try_from(record: FastqRecord<&'a [u8]>) -> Result<Self, Self::Error> {
        Ok(Self {
            buf: record.buf.try_borrow_str()?,
            ranges: record.ranges,
        })
    }
}

// StringBytes preserves the complete byte sequence when converting to text.
// Concrete sources avoid overlap with existing From conversions.
macro_rules! impl_record_try_into_string {
    ($($source:ty),+ $(,)?) => {
        $(
            impl TryFrom<FastqRecord<$source>> for FastqRecord<String> {
                type Error = FastqRecordError;

                /// Consume the record and convert its complete backing bytes to UTF-8 text.
                /// Borrowed slices are copied; vectors and boxed slices reuse their storage.
                /// Byte buffers reuse eligible unique storage and otherwise copy their bytes.
                /// Field ranges remain byte ranges and need not lie on character boundaries.
                ///
                /// # Errors
                /// Returns [`FastqRecordError::InvalidUtf8`] if any backing bytes are
                /// invalid UTF-8, including bytes outside the field ranges.
                /// The error retains the rejected bytes.
                fn try_from(record: FastqRecord<$source>) -> Result<Self, Self::Error> {
                    Ok(Self {
                        buf: record.buf.try_into_string()?,
                        ranges: record.ranges,
                    })
                }
            }
        )+
    };
}

impl_record_try_into_string!(&[u8], Vec<u8>, Box<[u8]>, Bytes, BytesMut);

/// A FASTQ record could not be validated or converted to text.
#[derive(Debug, Error)]
pub enum FastqRecordError {
    /// A field range is reversed or extends beyond the backing buffer.
    #[error("FASTQ {field} range is reversed or outside the backing buffer")]
    InvalidRange { field: &'static str },
    /// A field violates the FASTQ field constraints.
    #[error(transparent)]
    InvalidFields(#[from] FastqFieldError),
    /// The backing storage is invalid UTF-8; retains the rejected bytes.
    #[error(transparent)]
    InvalidUtf8(#[from] FromUtf8Error),
    /// The borrowed backing storage is invalid UTF-8.
    #[error(transparent)]
    InvalidUtf8Slice(#[from] Utf8Error),
}

// Keep the targets concrete: their vector conversions preserve every byte and
// therefore the newly constructed ranges. An arbitrary From<Vec<u8>> need not.
macro_rules! impl_record_from_field {
    ($($target:ty),+ $(,)?) => {
        $(
            impl<T: StringBytes> From<FastqField<T>> for FastqRecord<$target> {
                /// Allocate storage and concatenate the field bytes in FASTQ order.
                /// No format markers, spaces or line endings are inserted.
                fn from(fields: FastqField<T>) -> Self {
                    FastqRecord::<Vec<u8>>::from(fields).into()
                }
            }
        )+
    };
}

impl_record_from_field!(Box<[u8]>, Bytes, BytesMut);

impl<T: StringBytes> From<FastqField<T>> for FastqRecord<Vec<u8>> {
    /// Allocate a vector and concatenate the field bytes in FASTQ order.
    /// No format markers, spaces or line endings are inserted.
    fn from(fields: FastqField<T>) -> Self {
        let parts = fields.byte_view().into_parts();
        let capacity = parts.id.len()
            + parts.description.map_or(0, |desc| desc.len())
            + parts.sequence.len()
            + parts.separator.len()
            + parts.quality.len();
        let mut buf = Vec::with_capacity(capacity);
        let id = append(&mut buf, parts.id);
        let description = parts.description.map(|desc| append(&mut buf, desc));
        let sequence = append(&mut buf, parts.sequence);
        let separator = append(&mut buf, parts.separator);
        let quality = append(&mut buf, parts.quality);
        Self {
            buf,
            ranges: FastqParts {
                id,
                description,
                sequence,
                separator,
                quality,
            },
        }
    }
}

fn append(buf: &mut Vec<u8>, field: &[u8]) -> Range<usize> {
    let start = buf.len();
    buf.extend_from_slice(field);
    start..buf.len()
}

impl<T: StringBytes> TryFrom<(T, FastqParts<Range<usize>>)> for FastqRecord<T> {
    type Error = FastqRecordError;

    /// Validate the selected fields and retain the buffer without copying its bytes.
    /// Field ranges may overlap or appear in any order. Bytes outside the ranges
    /// are ignored; no `@` or `+` markers, field delimiters or line endings are required.
    /// The separator range selects only the content after the line's leading `+`.
    ///
    /// # Errors
    /// Returns [`FastqRecordError::InvalidRange`] for reversed or out-of-bounds
    /// ranges. Selected field contents are checked with [`FastqField::try_from`].
    fn try_from((buf, ranges): (T, FastqParts<Range<usize>>)) -> Result<Self, Self::Error> {
        let bytes = buf.as_bytes();
        let get = |field, range: &Range<usize>| {
            bytes
                .get(range.clone())
                .ok_or(FastqRecordError::InvalidRange { field })
        };
        let parts = FastqParts {
            id: get("identifier", &ranges.id)?,
            description: ranges
                .description
                .as_ref()
                .map(|range| get("description", range))
                .transpose()?,
            sequence: get("sequence", &ranges.sequence)?,
            separator: get("separator", &ranges.separator)?,
            quality: get("quality", &ranges.quality)?,
        };
        FastqField::try_from(parts)?;
        Ok(Self { buf, ranges })
    }
}

impl<T, U> TryFrom<FastqParts<T>> for FastqRecord<U>
where
    FastqField<T>: TryFrom<FastqParts<T>, Error = FastqFieldError>,
    FastqRecord<U>: From<FastqField<T>>,
{
    type Error = FastqFieldError;

    /// Check the parts, allocate storage and concatenate the field bytes in FASTQ order.
    /// No format markers, spaces or line endings are inserted.
    /// The target storage can be [`Vec<u8>`], `Box<[u8]>`, [`Bytes`] or [`BytesMut`].
    ///
    /// ```
    /// use biorill_fastx::{FastqParts, FastqRecord};
    ///
    /// let record = FastqRecord::<Vec<u8>>::try_from(FastqParts::new(
    ///     "read1",
    ///     Some("note"),
    ///     "AC",
    ///     "",
    ///     "12",
    /// ))?;
    /// assert_eq!(FastqRecord::<&[u8]>::from(&record).into_parts().0, b"read1noteAC12");
    /// let mut encoded = Vec::new();
    /// record.encode_into_vec(&mut encoded);
    /// assert_eq!(encoded, b"@read1 note\nAC\n+\n12\n");
    /// # Ok::<(), biorill_fastx::FastqFieldError>(())
    /// ```
    ///
    /// # Errors
    /// Returns the same errors as [`FastqField::try_from`].
    fn try_from(parts: FastqParts<T>) -> Result<Self, Self::Error> {
        FastqField::<T>::try_from(parts).map(Self::from)
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::ops::Range;
    use std::rc::Rc;
    use std::sync::Arc;

    use bytes::{Bytes, BytesMut};
    use stringbytes::StringBytes;

    use super::{FastqRecord, FastqRecordError};
    use crate::tests::valid_parts;
    use crate::{FastqField, FastqFieldError, FastqParts, LengthMismatch};

    fn parts() -> FastqParts<&'static [u8]> {
        FastqParts {
            id: b"id",
            description: Some(b"desc"),
            sequence: b"AC",
            separator: b"id",
            quality: b"12",
        }
    }

    #[test]
    fn owned_storage_variants_preserve_fields_and_bytes_conversion_keeps_the_buffer() {
        let fields = FastqField::try_from(parts()).unwrap();
        let mutable = FastqRecord::<BytesMut>::from(fields.clone());
        let pointer = FastqRecord::<&[u8]>::from(&mutable).into_parts().0.as_ptr();
        let (_, ranges) = FastqRecord::<&[u8]>::from(&mutable).into_parts();
        let frozen = FastqRecord::<Bytes>::from(mutable);
        assert_eq!(
            FastqRecord::<&[u8]>::from(&frozen).into_parts().0.as_ptr(),
            pointer
        );
        assert_eq!(frozen.parts(), parts());
        let (data, frozen_ranges) = frozen.into_parts();
        assert_eq!(frozen_ranges, ranges);

        let vector = FastqRecord::<Vec<u8>>::from(fields.clone());
        assert_eq!(vector.parts(), parts());
        assert_eq!(
            FastqRecord::<&[u8]>::from(&vector).into_parts().0,
            data.as_ref()
        );
        let (vector, vector_ranges) = vector.into_parts();
        assert_eq!(vector_ranges, ranges);
        assert_eq!(vector, b"iddescACid12");

        let boxed = FastqRecord::<Box<[u8]>>::from(fields.clone());
        assert_eq!(boxed.parts(), parts());
        let (boxed, boxed_ranges) = boxed.into_parts();
        assert_eq!(&*boxed, vector.as_slice());
        assert_eq!(boxed_ranges, ranges);

        let bytes = FastqRecord::<Bytes>::from(fields);
        assert_eq!(bytes.parts(), parts());
        assert_eq!(
            FastqRecord::<&[u8]>::from(&bytes).into_parts().0,
            vector.as_slice()
        );
    }

    #[test]
    fn checked_buffer_conversion_reuses_data_and_rejects_invalid_ranges() {
        let (data, ranges) =
            FastqRecord::<BytesMut>::from(FastqField::try_from(parts()).unwrap()).into_parts();
        for field in ["id", "description", "sequence", "separator", "quality"] {
            for range in [usize::MAX..usize::MAX, 0..data.len() + 1, 2..1] {
                let mut invalid = ranges.clone();
                match field {
                    "id" => invalid.id = range,
                    "description" => invalid.description = Some(range),
                    "sequence" => invalid.sequence = range,
                    "separator" => invalid.separator = range,
                    "quality" => invalid.quality = range,
                    _ => unreachable!(),
                }
                assert!(matches!(
                    FastqRecord::<BytesMut>::try_from((data.clone(), invalid)),
                    Err(FastqRecordError::InvalidRange { .. })
                ));
            }
        }
        let mut without_description = ranges.clone();
        without_description.description = None;
        let checked =
            FastqRecord::<BytesMut>::try_from((data.clone(), without_description)).unwrap();
        assert_eq!(checked.description(), None);
        assert_eq!(
            {
                let mut encoded = Vec::new();
                checked.encode_into_vec(&mut encoded);
                encoded
            },
            b"@id\nAC\n+id\n12\n"
        );
        let mut trailing = data.clone();
        trailing.extend_from_slice(b"extra");
        let checked = FastqRecord::<BytesMut>::try_from((trailing, ranges.clone())).unwrap();
        assert!(FastqRecord::<&[u8]>::from(&checked)
            .into_parts()
            .0
            .ends_with(b"extra"));
        assert_eq!(
            {
                let mut encoded = Vec::new();
                checked.encode_into_vec(&mut encoded);
                encoded
            },
            b"@id desc\nAC\n+id\n12\n"
        );

        let pointer = data.as_ptr();
        let checked = FastqRecord::<BytesMut>::try_from((data, ranges.clone())).unwrap();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&checked).into_parts().0.as_ptr(),
            pointer
        );
        assert_eq!(checked.parts(), parts());
        assert_eq!(checked.into_parts().1, ranges);
    }

    #[test]
    fn checked_buffer_conversion_ignores_unused_bytes_and_checks_field_contents() {
        let data = BytesMut::from(b"@id desc\nAC\n+id\n12\n".as_slice());
        let ranges = FastqParts::new(1..3, Some(4..8), 9..11, 13..15, 16..18);
        for offset in [
            0,
            ranges.id.end,
            ranges.description.as_ref().unwrap().end,
            ranges.sequence.end,
            ranges.separator.start - 1,
            ranges.separator.end,
            ranges.quality.end,
        ] {
            let mut invalid = data.clone();
            invalid[offset] = b'?';
            let checked = FastqRecord::<BytesMut>::try_from((invalid, ranges.clone())).unwrap();
            assert_eq!(checked.parts(), parts());
            assert_eq!(
                {
                    let mut encoded = Vec::new();
                    checked.encode_into_vec(&mut encoded);
                    encoded
                },
                data.as_ref()
            );
        }
        for (offset, value) in [
            (ranges.id.start, b' '),
            (ranges.description.as_ref().unwrap().start, b'\n'),
            (ranges.sequence.start, b'\n'),
            (ranges.separator.start, b'\n'),
            (ranges.quality.start, b'\n'),
        ] {
            let mut invalid = data.clone();
            invalid[offset] = value;
            assert!(matches!(
                FastqRecord::<BytesMut>::try_from((invalid, ranges.clone())),
                Err(FastqRecordError::InvalidFields(_))
            ));
        }
        let mut short_quality = data.to_vec();
        short_quality.remove(ranges.quality.start);
        let mut short_ranges = ranges;
        short_ranges.quality.end -= 1;
        assert!(matches!(
            FastqRecord::<BytesMut>::try_from((
                BytesMut::from(short_quality.as_slice()),
                short_ranges
            )),
            Err(FastqRecordError::InvalidFields(
                FastqFieldError::LengthMismatch(LengthMismatch { seq: 2, qual: 1 })
            ))
        ));
    }

    #[test]
    fn parts_borrow_the_complete_record_and_edits_leave_it_unchanged() {
        let original = FastqRecord::<Bytes>::try_from(parts()).unwrap();
        let mut parts = original.parts();
        for (field, start) in [
            (parts.id, 0),
            (parts.description.unwrap(), 2),
            (parts.sequence, 6),
            (parts.separator, 8),
            (parts.quality, 10),
        ] {
            assert_eq!(
                field.as_ptr(),
                FastqRecord::<&[u8]>::from(&original).into_parts().0[start..].as_ptr()
            );
        }
        let borrowed = FastqRecord::<&[u8]>::from(&original);
        assert_eq!(
            FastqRecord::<&[u8]>::from(&borrowed)
                .into_parts()
                .0
                .as_ptr(),
            FastqRecord::<&[u8]>::from(&original)
                .into_parts()
                .0
                .as_ptr()
        );
        parts.id = b"new";
        parts.description = None;
        parts.sequence = b"C";
        parts.separator = b"new";
        parts.quality = b"2";
        let updated = FastqRecord::<Bytes>::try_from(parts).unwrap();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&updated).into_parts().0,
            b"newCnew2"
        );
        assert_eq!(
            FastqRecord::<&[u8]>::from(&original).into_parts().0,
            b"iddescACid12"
        );
    }

    #[test]
    fn invalid_parts_fail_without_changing_the_original_record() {
        fn check<U: StringBytes + AsRef<[u8]>>()
        where
            for<'a> FastqRecord<U>: From<FastqField<&'a [u8]>>,
        {
            let original = FastqRecord::<U>::try_from(parts()).unwrap();
            let mut parts = original.parts();
            parts.quality = b"1";
            assert!(matches!(
                FastqRecord::<U>::try_from(parts),
                Err(FastqFieldError::LengthMismatch(LengthMismatch {
                    seq: 2,
                    qual: 1
                }))
            ));
            let mut parts = original.parts();
            parts.id = b"id extra";
            assert!(matches!(
                FastqRecord::<U>::try_from(parts),
                Err(FastqFieldError::InvalidId)
            ));
            let mut parts = original.parts();
            parts.separator = b"id\nextra";
            assert!(matches!(
                FastqRecord::<U>::try_from(parts),
                Err(FastqFieldError::LineFeed { field: "separator" })
            ));
            let mut parts = original.parts();
            parts.description = Some(b"extra\nline");
            assert!(matches!(
                FastqRecord::<U>::try_from(parts),
                Err(FastqFieldError::LineFeed {
                    field: "description"
                })
            ));
            assert_eq!(
                FastqRecord::<&[u8]>::from(&original).into_parts().0,
                b"iddescACid12"
            );
        }

        check::<Vec<u8>>();
        check::<Box<[u8]>>();
        check::<Bytes>();
        check::<BytesMut>();
    }

    #[test]
    fn constructing_records_reindexes_growing_shrinking_and_empty_parts() {
        fn check<U: StringBytes + AsRef<[u8]>>()
        where
            for<'a> FastqRecord<U>: From<FastqField<&'a [u8]>>,
        {
            let original = FastqRecord::<U>::try_from(parts()).unwrap();
            for description in [
                Some(b"long description".as_slice()),
                Some(b"x".as_slice()),
                Some(b"".as_slice()),
                None,
            ] {
                for (sequence, quality) in [
                    (b"ACGTA".as_slice(), b"12345".as_slice()),
                    (b"T".as_slice(), b"5".as_slice()),
                    (b"".as_slice(), b"".as_slice()),
                ] {
                    let mut parts = original.parts();
                    parts.description = description;
                    parts.sequence = sequence;
                    parts.quality = quality;
                    let updated = FastqRecord::<U>::try_from(parts.clone()).unwrap();
                    assert_eq!(updated.id(), b"id");
                    assert_eq!(updated.description(), description);
                    assert_eq!(updated.sequence(), sequence);
                    assert_eq!(updated.separator(), b"id");
                    assert_eq!(updated.quality(), quality);
                    assert_eq!(updated.parts(), parts);
                    assert_eq!(
                        FastqRecord::<&[u8]>::from(&updated).into_parts().0,
                        [
                            b"id".as_slice(),
                            description.unwrap_or_default(),
                            sequence,
                            b"id",
                            quality,
                        ]
                        .concat()
                    );
                }
            }
        }

        check::<Vec<u8>>();
        check::<Box<[u8]>>();
        check::<Bytes>();
        check::<BytesMut>();
    }

    #[test]
    fn owned_parts_build_a_record_and_preserve_arbitrary_sequence_bytes() {
        let record = FastqRecord::<Bytes>::try_from(FastqParts {
            id: Vec::new(),
            description: None,
            sequence: vec![b'A', 255],
            separator: b"".to_vec(),
            quality: b"12".to_vec(),
        })
        .unwrap();
        assert_eq!(
            FastqRecord::<&[u8]>::from(&record).into_parts().0,
            b"A\xff12"
        );
        assert_eq!(
            {
                let mut encoded = Vec::new();
                record.encode_into_vec(&mut encoded);
                encoded
            },
            b"@\nA\xff\n+\n12\n"
        );
    }

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
        assert!(fields.set_sequence_and_quality(b"AC", b"1").is_err());
        assert_eq!(fields, original);

        fields.set_id(b"read1").unwrap();
        fields.set_description(Some(b"note")).unwrap();
        fields.set_separator(b"").unwrap();
        fields.set_sequence_and_quality(b"AC", b"12").unwrap();
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
            let text =
                <FastqRecord<String> as TryFrom<(String, FastqParts<Range<usize>>)>>::try_from((
                    data,
                    ranges(),
                ))
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
            assert_eq!(record.parts(), valid_parts());
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
    fn record_text_views_validate_selected_fields_and_reject_split_characters() {
        let record = FastqRecord::try_from((
            b"\xffread1AC12\xff".as_slice(),
            FastqParts::new(1..6, None, 6..8, 8..8, 8..10),
        ))
        .unwrap();
        let fields = record.fields();
        assert_eq!(fields, FastqField::from(&record));
        let text = fields.try_str_view().unwrap();
        assert_eq!(text.to_string(), "@read1\nAC\n+\n12\n");
        assert_eq!(text.id().as_ptr(), record.id().as_ptr());
        assert_eq!(text.sequence().as_ptr(), record.sequence().as_ptr());
        assert_eq!(text.quality().as_ptr(), record.quality().as_ptr());

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
            fields.try_str_view(),
            Err(FastqFieldError::InvalidUtf8Slice {
                field: "sequence",
                ..
            })
        ));
    }
}
