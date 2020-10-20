use std::path::Path;
use std::error::Error;
use std::fs::File;

use structs1::{SegRec, LinkRec, HeaderRec};
use prepack1::{Gfa1Prepack, RecordType};

// use ::gfa_graph::parser::prepack1::RecordType;

pub mod prepack1;
pub mod opt_fields;
pub mod structs1;

pub fn parse_gfa_v1(gfa_file: &Path) -> Result<Gfa1Prepack, Box<dyn Error>> {
    let mut prepack = Gfa1Prepack::new();

    let file = File::open(gfa_file)?;
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .comment(Some(b'#'))
        .from_reader(file);

    let mut raw_record = csv::ByteRecord::new();
    while rdr.read_byte_record(&mut raw_record)? {
        if let Some(rec_type) = RecordType::from_raw(&raw_record[0]) {
            match rec_type {
                RecordType::Comment => {},
                RecordType::Header => {prepack.update_header(HeaderRec::from_raw(&raw_record)?); },
                RecordType::Sequence => { prepack.add_segment(SegRec::from_raw(&raw_record)?); },
                RecordType::Link => { prepack.add_link(LinkRec::from_raw(&raw_record)?); }
                _ => {println!("Support of this record type would be added in future.")}
            };
        } else {
            println!("Unsupported record type. We ignore it.")
        }
    }

    Ok(prepack)
}

//pub fn parse_gfa_v2(gfa_file: &Path) -> Result<Gfa1Prepack, Box<dyn Error>> {
//}