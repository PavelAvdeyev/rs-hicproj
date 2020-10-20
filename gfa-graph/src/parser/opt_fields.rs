use lazy_static::lazy_static;
use regex::bytes::Regex;
use std::str;
use ascii::AsciiString;

pub type OptFieldTag = AsciiString;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum OptFieldVal {
    A(u8),
    Int(i64),
    Float(f64),
    H(Vec<u32>),
    Z(AsciiString),
    // J(AsciiString),
    // BInt(Vec<i64>),
    // BFloat(Vec<f32>),
}

impl OptFieldVal {
    pub fn parse(input: &[u8]) -> Option<OptFieldVal> {
        lazy_static! {
            static ref RE_TAG: Regex = Regex::new(r"(?-u)[A-Za-z][A-Za-z0-9]").unwrap();
            static ref RE_CHAR: Regex = Regex::new(r"(?-u)[!-~]").unwrap();
            static ref RE_INT: Regex = Regex::new(r"(?-u)[-+]?[0-9]+").unwrap();
            static ref RE_FLOAT: Regex = Regex::new(r"(?-u)[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?")
                    .unwrap();
            static ref RE_STRING: Regex = Regex::new(r"(?-u)[ !-~]+").unwrap();
            static ref RE_BYTES: Regex = Regex::new(r"(?-u)[0-9A-F]+").unwrap();
        }

        let o_type = input.get(0)?;
        if !b"AifZJHB".contains(&o_type) {
            return None;
        }

        let o_contents = input.get(2..)?;
        match o_type {
            // char
            b'A' => RE_CHAR.find(o_contents).map(|s| s.as_bytes()[0]).map(OptFieldVal::A),
            // int
            b'i' => RE_INT
                .find(o_contents)
                .and_then(|s| str::from_utf8(s.as_bytes()).ok())
                .and_then(|s| s.parse().ok())
                .map(OptFieldVal::Int),
            // float
            b'f' => RE_FLOAT
                .find(o_contents)
                .and_then(|s| str::from_utf8(s.as_bytes()).ok())
                .and_then(|s| s.parse().ok())
                .map(OptFieldVal::Float),
            // bytearray
            b'H' => RE_BYTES
                .find(o_contents)
                .and_then(|s| str::from_utf8(s.as_bytes()).ok())
                .map(|s| s.chars().filter_map(|c| c.to_digit(16)))
                .map(|s| OptFieldVal::H(s.collect())),
            b'Z' => RE_STRING
                .find(o_contents)
                .map(|s| OptFieldVal::Z(AsciiString::from_ascii(s.as_bytes())
                    .expect("Problem with parsing Z tag. Non-ascii character is present"))),
            _ => panic!(
                "Tried to parse optional field with unknown type '{}'",
                o_type,
            ),
        }
    }
}

impl std::fmt::Display for OptFieldVal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OptFieldVal::A(x) => write!(f, "A:{}", char::from(*x)),
            OptFieldVal::Int(x) => write!(f, "i:{}", *x),
            OptFieldVal::Float(x) => write!(f, "f:{}", *x),
            OptFieldVal::H(x) => {
                write!(f, "H:")?;
                for a in x {
                    write!(f, "{:x}", *a)?
                }
                Ok(())
            }
            OptFieldVal::Z(x) => write!(f, "Z:{}", x),
        }
    }
}

/// Parses an optional field from a bytestring in the format <TAG>:<TYPE>:<VALUE>
pub fn parse_opt_field(input: &[u8]) -> Option<(OptFieldTag, OptFieldVal)> {
    let o_tag = match convert_to_tag(input.get(0..=1)?) {
        Some(x) => x,
        _ => return None,
    };

    let o_val = match OptFieldVal::parse(input.get(3..)?) {
        Some(x) => x,
        _ => return None,
    };

    Some((o_tag, o_val))
}

fn convert_to_tag(t: &[u8]) -> Option<OptFieldTag> {
    if t.len() != 2 || !t[0].is_ascii_alphabetic() || !t[1].is_ascii_alphanumeric() {
        return None;
    }
    OptFieldTag::from_ascii(t.clone()).ok()
}

// impl std::fmt::Display for OptField {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         use OptFieldVal::*;
//
//         write!(f, "{}{}:", char::from(self.tag[0]), char::from(self.tag[1]))?;
//
//         match &self.value {
//             A(x) => write!(f, "A:{}", char::from(*x)),
//             Int(x) => write!(f, "i:{}", x),
//             Float(x) => write!(f, "f:{}", x),
//             H(x) => {
//                 write!(f, "H:")?;
//                 for a in x {
//                     write!(f, "{:x}", a)?
//                 }
//                 Ok(())
//             }
//             // Z(x) => write!(f, "Z:{}", x),
//             // J(x) => write!(f, "J:{}", x),
//             // BInt(x) => {
//             //     write!(f, "B:I{}", x[0])?;
//             //     for a in x[1..].iter() {
//             //         write!(f, ",{}", a)?
//             //     }
//             //     Ok(())
//             // }
//             // BFloat(x) => {
//             //     write!(f, "B:F{}", x[0])?;
//             //     for a in x[1..].iter() {
//             //         write!(f, ",{}", a)?
//             //     }
//             //     Ok(())
//             // }
//         }
//     }
// }

