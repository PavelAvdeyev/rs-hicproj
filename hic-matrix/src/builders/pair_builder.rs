use ahash::AHashMap;
use std::path::{Path, PathBuf};
use ascii::{AsciiString, AsAsciiStr, AsAsciiStrError, AsciiStr};
use ndarray::{Array1, ArrayView1};
use std::iter::FromIterator;
use std::error::Error;
use std::fs::File;
use serde::Deserialize;

use super::res_grp_builder::ResGrpBuilder;
use super::super::reader::PixelT;

pub struct PairsBuilder {
    rsltn: u32,
    n_bins: usize,
    name2order: AHashMap<AsciiString, usize>,
    tig_order: Array1<AsciiString>,
    tig_lengths: Array1<u64>,
    bin_table: (Array1<u32>, Array1<u64>, Array1<u64>),
    tig_offsets: Array1<u32>,
    pairs_file: PathBuf,
}

#[derive(Debug, Deserialize)]
struct PairRecord<'a> {
    read_name: &'a str,
    tig1: &'a str,
    pos1: u64,
    tig2: &'a str,
    pos2: u64,
    strand1: char,
    strand2: char,
}

impl ResGrpBuilder for PairsBuilder {
    fn get_resolution(&self) -> u32 {
        self.rsltn
    }

    fn get_tig_offsets_view(&self) -> ArrayView1<u32> {
        self.tig_offsets.view()
    }

    fn get_bin_table(&self) -> (ArrayView1<u32>, ArrayView1<u64>, ArrayView1<u64>) {
        (self.bin_table.0.view(), self.bin_table.1.view(), self.bin_table.2.view())
    }

    fn get_bin_offsets(&self, pixels: &[PixelT]) -> Array1<u32> {
        PairsBuilder::build_bin_offsets_from_pixels(self.n_bins, pixels)
    }

    fn get_pixels(&self) -> Result<Vec<PixelT>, Box<dyn Error>> {
        let mut pixels:AHashMap<(u32, u32), u32> = AHashMap::default();
        let file = File::open(self.pairs_file.as_path())?;

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .has_headers(false)
            .from_reader(file);
        let mut raw_record = csv::ByteRecord::new();
        let mut total: u32 = 0;

        while rdr.read_byte_record(&mut raw_record)? {
            total += 1;
            let record: PairRecord = raw_record.deserialize(None)?;

            if let Some(bin_rec) = self.pair_to_bin_rec(&record)? {
                let count = pixels.entry(bin_rec).or_insert(0);
                *count += 1;
            } else {
                println!("There is a problem with pair record.");
            }

            if total % 1000000 == 0 {
                println!("{} hic pairs were converted to bins", total);
            }
        }

        let mut pixels: Vec<PixelT> = pixels.into_iter().map(|x| (x.0.0, x.0.1, x.1)).collect();
        pixels.sort_by_key(|rec| { (rec.0, rec.1) });

        Ok(pixels)
    }

}

impl PairsBuilder {
    pub fn new(pairs_file: &Path, ord_tig_lengths: &[(AsciiString, u64)], rsltn: u32) -> PairsBuilder {
        let tig_lengths: Array1<u64> = Array1::from_iter(ord_tig_lengths.iter().map(|x| x.1));
        let tig_offsets = PairsBuilder::build_tig_offsets(rsltn,tig_lengths.view());
        let n_bins = if !tig_offsets.is_empty() {tig_offsets[tig_offsets.len() - 1] as usize} else {0};
        let bin_table = PairsBuilder::build_bin_table_from_lengths(n_bins, rsltn as u64, tig_lengths.view());

        PairsBuilder {
            rsltn,
            n_bins,
            name2order: ord_tig_lengths.iter().enumerate()
                .map(|(i, x)| (x.0.clone(), i) )
                .collect(),
            tig_order: Array1::from_iter(ord_tig_lengths.iter()
                .map(|x| x.0.clone())),
            tig_lengths,
            bin_table,
            tig_offsets,
            pairs_file: PathBuf::from(pairs_file),
        }
    }

    pub fn tig_names_view(&self) -> ArrayView1<AsciiString> {
        self.tig_order.view()
    }

    pub fn tig_lengths_view(&self) -> ArrayView1<u64> {
        self.tig_lengths.view()
    }

    fn pair_to_bin_rec(&self, record: &PairRecord) -> Result<Option<(u32, u32)>, AsAsciiStrError> {
        let get_anchor = |anchor: u64, tig: &AsciiStr| -> Option<u64> {
            self.get_tig_length_by_name(tig).map(|max_len| max_len.min(anchor))
        };

        let tig1 = AsciiString::from(record.tig1.as_ascii_str()?);
        let tig2 = AsciiString::from(record.tig2.as_ascii_str()?);
        let rsltn = self.rsltn as u64;

        let anchors = get_anchor(record.pos1, &tig1).zip(get_anchor(record.pos2, &tig2));
        let offsets = self.get_tig_offset_by_name(&tig1).zip(self.get_tig_offset_by_name(&tig2));
        let bin_ids = offsets.zip(anchors)
            .map(|(offset, anchor)|
                (offset.0 + (anchor.0 / rsltn) as u32, offset.1 + (anchor.1 / rsltn) as u32)
            )
            .map(|(bin1, bin2)| if bin1 <= bin2 { (bin1, bin2) } else { (bin2, bin1) });
        Ok(bin_ids)
    }

    fn get_tig_length_by_name(&self, tig: &AsciiStr) -> Option<u64> {
        self.name2order.get(tig).map(|id| self.tig_lengths[*id])
    }

    fn get_tig_offset_by_name(&self, tig: &AsciiStr) -> Option<u32> {
        self.name2order.get(tig).map(|id| self.tig_offsets[*id])
    }
}

