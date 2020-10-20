use ascii::AsciiString;
use std::fmt;
use std::str::FromStr;

const FIELD_SEP: char = '\t';
// pub const COL_READID: usize = 0;
pub const COL_TIG1: usize = 1;
pub const COL_POS1: usize = 2;
pub const COL_TIG2: usize = 3;
pub const COL_POS2: usize = 4;
// pub const COL_STRAND1: usize = 5;
// pub const COL_STRAND2: usize = 6;

pub enum Strand {
    Forward,
    Reverse
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-")
        }
    }
}

pub struct PairRecord {
    pub qname: AsciiString,
    pub name1: AsciiString,
    pub pos1: i64,
    pub strand1: Strand,
    pub name2: AsciiString,
    pub pos2: i64,
    pub strand2: Strand
}

impl PairRecord {
    pub fn _new() -> PairRecord {
        PairRecord {
            qname: AsciiString::default(),
            name1: AsciiString::default(),
            pos1: -1,
            strand1: Strand::Forward,
            name2: AsciiString::default(),
            pos2: -1,
            strand2: Strand::Forward,
        }
    }

    pub fn from_bams(r1: &bam::Record, r2: &bam::Record) -> PairRecord {
        let (r1, r2) = get_ordered_alignments(r1, r2);
        PairRecord {
            qname: AsciiString::from_ascii(r1.name()).unwrap(),
            name1: AsciiString::from_ascii(r1.ref_id().to_string().as_bytes()).unwrap(),
            pos1: get_alignment_pos(r1),
            strand1: if r1.flag().is_reverse_strand() {Strand::Reverse} else {Strand::Forward},
            name2: AsciiString::from_ascii(r2.ref_id().to_string().as_bytes()).unwrap(),
            pos2: get_alignment_pos(r2),
            strand2: if r2.flag().is_reverse_strand() {Strand::Reverse} else {Strand::Forward},
        }
    }

    // pub fn get_ordered_coordinates(&self) -> Option<(i64, i64)> {
    //     if self.name1 == self.name2 {
    //         if self.pos1 < self.pos2 {
    //             Some((self.pos1, self.pos2))
    //         } else {
    //             Some((self.pos2, self.pos1))
    //         }
    //     } else { None }
    // }

    pub fn to_string(&self) -> AsciiString {
        AsciiString::from_str(
            format!("{0}{7}{1}{7}{2}{7}{3}{7}{4}{7}{5}{7}{6}", self.qname, self.name1, self.pos1, self.name2, self.pos2, self.strand1, self.strand2, FIELD_SEP).as_str()
        ).unwrap()
    }
}

pub fn get_alignment_pos(rec: &bam::Record) -> i64 {
    (rec.start() + (rec.calculate_end() - rec.start()) / 2) as i64
}

pub fn calc_dist(fa: &bam::Record, sa: &bam::Record) -> u64 {
    if fa.ref_id() == sa.ref_id() {
        (get_alignment_pos(fa) - get_alignment_pos(sa)).abs() as u64
    } else {
        u64::MAX
    }
}

pub fn is_opposite_pair(fa: &bam::Record, sa: &bam::Record) -> bool {
    let mut can_rescue = fa.flag().is_reverse_strand() != sa.flag().is_reverse_strand();
    if sa.flag().is_reverse_strand() {
        can_rescue &= get_alignment_pos(sa) > get_alignment_pos(fa)
    } else {
        can_rescue &= get_alignment_pos(sa) < get_alignment_pos(fa)
    }
    can_rescue
}

pub fn get_ordered_alignments<'a>(fa: &'a bam::Record, sa: &'a bam::Record) -> (&'a bam::Record, &'a bam::Record) {
    if fa.ref_id() == sa.ref_id() {
        if get_alignment_pos(fa) < get_alignment_pos(sa) { (fa, sa) } else { (sa, fa) }
    } else if fa.ref_id() < sa.ref_id() {
        (fa, sa)
    } else {
        (sa, fa)
    }
}



