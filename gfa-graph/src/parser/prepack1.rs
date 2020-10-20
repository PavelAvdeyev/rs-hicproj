use super::structs1::{SegRec, LinkRec, HeaderRec};

pub enum RecordType {
    Comment,
    Header,
    Sequence,
    Link,
    Containment,
    Path
}

impl RecordType {
    pub fn from_raw(s: &[u8]) -> Option<RecordType> {
        match *s {
            [b'#'] => Some(RecordType::Comment),
            [b'H'] => Some(RecordType::Header),
            [b'S'] => Some(RecordType::Sequence),
            [b'L'] => Some(RecordType::Link),
            [b'C'] => Some(RecordType::Containment),
            [b'P'] => Some(RecordType::Path),
            _ => None
        }
    }
}

// #[derive(Default, Debug, Clone)]
pub struct Gfa1Prepack {
    header: HeaderRec,
    sequences: Vec<SegRec>,
    links: Vec<LinkRec>,
}

impl Gfa1Prepack {
    pub fn new() -> Gfa1Prepack {
        Gfa1Prepack {
            header: Default::default(),
            sequences: Vec::new(),
            links: Vec::new()
        }
    }

    pub fn from(sequences: Vec<SegRec>, links: Vec<LinkRec>) -> Gfa1Prepack {
        Gfa1Prepack {
            header: Default::default(),
            sequences,
            links
        }
    }

    pub fn update_header(&mut self, rec: HeaderRec) {
        self.header = rec;
    }

    pub fn add_segment(&mut self, rec: SegRec) {
        self.sequences.push(rec);
    }

    pub fn add_link(&mut self, rec: LinkRec) {
        self.links.push(rec);
    }

    pub fn seq_recs_iter<'a>(&'a self) -> impl Iterator<Item = &SegRec> + 'a {
        self.sequences.iter()
    }

    pub fn link_recs_iter<'a>(&'a self) -> impl Iterator<Item = &LinkRec> + 'a{
        self.links.iter()
    }
}

