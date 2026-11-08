mod field;
mod parts;
mod record;

pub use self::field::{FastqField, FastqFieldError, LengthMismatch};
pub use self::parts::FastqParts;
pub use self::record::{FastqRecord, FastqRecordError};
