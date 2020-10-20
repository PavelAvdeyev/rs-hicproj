use ascii::AsciiString;
use std::str::FromStr;
use std::ops::Add;
use std::fmt;

#[derive(Debug, Clone)]
pub enum Orientation {
    Forward,
    Reverse,
}


impl Orientation {
    pub fn from_raw(o: &[u8]) -> Option<Orientation> {
        match o {
            [b'+'] => Some(Orientation::Forward),
            [b'-'] => Some(Orientation::Reverse),
            _ => None
        }
    }

    pub fn inverse(o: &Orientation) -> Orientation {
        match o {
            Orientation::Forward => Orientation::Reverse,
            Orientation::Reverse => Orientation::Forward,
        }
    }

    pub fn to_node_name(s: AsciiString, o: Orientation) -> AsciiString {
        match o {
            Orientation::Forward => s.add(&AsciiString::from_str("+").unwrap()),
            Orientation::Reverse => s.add(&AsciiString::from_str("-").unwrap()),
        }
    }
}

impl fmt::Display for Orientation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Orientation::Forward => write!(f, "+"),
            Orientation::Reverse => write!(f, "-")
        }
    }
}
