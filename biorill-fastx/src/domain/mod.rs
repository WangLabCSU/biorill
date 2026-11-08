mod pair;
mod parser;
mod vo;

pub use self::pair::FastqPairError;
pub use self::parser::ParseError;
pub(crate) use self::parser::{FastqParser, ParsingStep, ReadingState};
pub use self::vo::{
    FastqField, FastqFieldError, FastqParts, FastqRecord, FastqRecordError, LengthMismatch,
};

#[cfg(test)]
pub(crate) use self::parser::RecordState;
