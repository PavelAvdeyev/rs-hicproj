use std::path::Path;
use std::fs::File;
use std::collections::VecDeque;
use std::io::{BufWriter, Write};

use log::info;
use serde::Deserialize;

pub fn deduplicate_pairs(inp_file: &Path, out_file: &Path, max_mismatch: i64) {
    // Find and remove PCR/optical duplicates.
    // Find PCR duplicates in an upper-triangular flipped sorted pairs file.
    // Allow for a +/-N bp mismatch at each side of duplicated molecules.
    let input = File::open(inp_file).unwrap();
    let output= File::create(out_file).unwrap();

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(input);
    let mut wrtr = BufWriter::new(output);

    let mut total: u32 = 0;
    let mut raw_record = csv::ByteRecord::new();
    let mut cur_records = VecDeque::<(Record, bool)>::new();
    while rdr.read_byte_record(&mut raw_record).unwrap() {
        total += 1;

        let low = cur_records.pop_front();
        let high: Record = raw_record.deserialize(None).unwrap();

        match low {
            None => cur_records.push_back((high, false)),
            Some((low, l_rm)) => {
                if l_rm {
                    cur_records.push_back((high, false));
                    cur_records = update_records_wrt_first(cur_records);
                } else if low.name1 != high.name1 || low.name2 != high.name2
                    || high.pos1 - low.pos1 > max_mismatch || high.pos1 - low.pos1 < 0
                    || high.pos2 - low.pos2 > max_mismatch || high.pos2 - low.pos2 < 0 {
                    // if we jumped too far, continue
                    save_record(&low, &mut wrtr);
                    cur_records.push_back((high, false));
                    cur_records = update_records_wrt_first(cur_records);
                } else if is_duplicated_copies(&low, &high, max_mismatch) {
                    cur_records.push_front((low, l_rm));
                } else {
                    cur_records.push_front((low, l_rm));
                    cur_records.push_back((high, false));
                }
            }
        }

        if total % 10000000 == 0 {
            info!("{} hic pairs were checked", total);
        }
    }
}

#[derive(Debug, Deserialize, Clone)]
struct Record {
    read_name: String,
    name1: String,
    pos1: i64,
    name2: String,
    pos2: i64,
    strand1: char,
    strand2: char,
}

fn is_duplicated_copies(rec1: &Record, rec2: &Record, max_mismatch: i64) -> bool {
    rec1.name1 == rec2.name1 && rec1.name2 == rec2.name2
        && rec1.strand1 == rec2.strand1 && rec1.strand2 == rec2.strand2
        && (rec1.pos1 - rec2.pos1).abs().max(rec1.pos2 - rec2.pos2) <= max_mismatch
}

fn save_record(rec: &Record, output: &mut BufWriter<File>) {
    writeln!(output, "{}", rec.read_name).expect("Problem with writing file");
}

fn update_records_wrt_first(mut recs: VecDeque<(Record, bool)>) -> VecDeque<(Record, bool)> {
    if recs.len() < 2 { return recs; }

    let f_rec = recs.front().unwrap().clone();
    if f_rec.1 { return recs; }

    for elem in recs.iter_mut().skip(1) {
        if is_duplicated_copies(&f_rec.0, &elem.0, 3) { elem.1 = true; }
    }
    recs
}
