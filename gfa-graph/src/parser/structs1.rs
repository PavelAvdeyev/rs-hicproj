use std::{fmt, error, str::FromStr};

use ascii::{AsciiString, AsAsciiStr};
use csv::ByteRecord;

use super::super::utils::Orientation;
use super::opt_fields::{self, OptFieldVal, OptFieldTag};

type OptionalFields = ahash::AHashMap<OptFieldTag, OptFieldVal>;
pub const OMITTED_SEQ_SYMBOL: u8 = b'*';


#[derive(Default, Debug, Clone)]
pub struct HeaderRec {
    optionals: OptionalFields
}

impl HeaderRec {
    pub fn from_raw(s: &ByteRecord) -> Result<HeaderRec, Box<dyn error::Error>> {
        Ok(HeaderRec{
            optionals: init_opt_fields(s, 1),
        })
    }

    pub fn get_version_number(&self) -> Option<&AsciiString>{
        get_tag_val("VN", &self.optionals).and_then(|x| match x {
            OptFieldVal::Z(val) => Some(val),
            _ => None
        })
    }

    pub fn get_tag_value(&self, name: &str) -> Option<&OptFieldVal> {
        get_tag_val(name, &self.optionals)
    }
}

#[derive(Clone)]
pub struct SegRec {
    pub name: AsciiString,
    pub seq: Option<AsciiString>,
    optionals: OptionalFields
}

impl SegRec {
    pub fn from_raw(s: &ByteRecord) -> Result<SegRec, Box<dyn error::Error>>{
        if s.len() < 3 {
            return Err(GFAParseError.into());
        }

        Ok(SegRec {
            name: AsciiString::from(s[1].as_ascii_str()?),
            seq: if s[2][0] != OMITTED_SEQ_SYMBOL { Some(AsciiString::from(s[2].as_ascii_str()?)) } else { None },
            optionals: init_opt_fields(s, 3),
        })
    }

    pub fn get_length(&self) -> Option<u64> {
        get_tag_val("LN", &self.optionals).and_then(|x| match x {
            OptFieldVal::Int(val) => Some(*val as u64),
            _ => None
        })
    }

    pub fn get_read_count(&self) -> Option<u64> {
        get_tag_val("RC", &self.optionals).and_then(|x| match x {
            OptFieldVal::Int(val) => Some(*val as u64),
            _ => None
        })
    }

    pub fn get_fragment_count(&self) -> Option<u64> {
        get_tag_val("FC", &self.optionals).and_then(|x| match x {
            OptFieldVal::Int(val) => Some(*val as u64),
            _ => None
        })
    }

    pub fn get_kmer_count(&self) -> Option<u64> {
        get_tag_val("KC", &self.optionals).and_then(|x| match x {
            OptFieldVal::Int(val) => Some(*val as u64),
            _ => None
        })
    }

    // SH	H	SHA-256 checksum of the sequence
    // pub fn get_sha_cheksum(&self) -> Option<&AsciiString> {
    //
    // }

    // UR	Z	URI or local file-system path of the sequence. If it does not start with a standard protocol (e.g. ftp), it is assumed to be a local path.
    // pub fn get_sequence_path(&self) -> Option<&AsciiString> {
    //
    // }

    pub fn get_tag_value(&self, name: &str) -> Option<&OptFieldVal> {
        get_tag_val(name, &self.optionals)
    }
}


#[derive(Debug, Clone)]
pub struct LinkRec {
    pub from_name: AsciiString,
    pub from_strand: Orientation,
    pub to_name: AsciiString,
    pub to_strand: Orientation,
    pub cigar: AsciiString,
    optionals: OptionalFields
}

impl LinkRec {
    pub fn from_raw(s: &ByteRecord) -> Result<LinkRec, Box<dyn error::Error>>{
        if s.len() < 6 {
            return Err(GFAParseError.into());
        }

        Ok(LinkRec {
            from_name: AsciiString::from(s[1].as_ascii_str()?),
            from_strand: Orientation::from_raw(&s[2]).ok_or(GFAParseError)?,
            to_name: AsciiString::from(s[3].as_ascii_str()?),
            to_strand:Orientation::from_raw(&s[4]).ok_or(GFAParseError)?,
            cigar: AsciiString::from(s[5].as_ascii_str()?),
            optionals: init_opt_fields(s, 6)
        })
    }

    pub fn inverse(l: &LinkRec) -> LinkRec {
        LinkRec {
            from_name: l.from_name.clone(),
            from_strand: Orientation::inverse(&l.from_strand),
            to_name: l.to_name.clone(),
            to_strand: Orientation::inverse(&l.to_strand),
            cigar: l.cigar.clone(),
            optionals: l.optionals.clone()
        }
    }
}

#[derive(Debug, Clone)]
pub struct ContainmentRec {
    pub container_name: AsciiString,
    pub container_orient: Orientation,
    pub contained_name: AsciiString,
    pub contained_orient: Orientation,
    pub pos: usize,
    pub overlap: AsciiString,
    optional: OptionalFields,
}

pub struct PathRec {
    pub path_name: AsciiString,
    //pub segment_names: AsciiString,
    //pub overlaps: Vec<BString>,
    optional: OptionalFields,
}


fn init_opt_fields(s: &ByteRecord, s_ind: usize) -> OptionalFields {
    let mut opts = OptionalFields::default();
    if s.len() > s_ind {
        for i in s_ind..s.len() {
            if let Some((tag, val)) = opt_fields::parse_opt_field(&s[i]) {
                opts.insert(tag, val);
            }
        }
    }
    opts
}

fn get_tag_val<'a>(name: &str, opt_fields: &'a OptionalFields) -> Option<&'a OptFieldVal> {
    AsciiString::from_str(name).ok().and_then(|key| opt_fields.get(&key))
}

#[derive(Debug, Clone)]
pub struct GFAParseError;

impl fmt::Display for GFAParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Problem with parsing file")
    }
}

impl error::Error for GFAParseError {}
