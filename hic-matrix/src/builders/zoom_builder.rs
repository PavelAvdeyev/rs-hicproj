use ndarray::{ArrayView1, Array1};
use std::error::Error;
use ahash::AHashMap;

use super::super::{res_group::ResGroup, reader::PixelT};
use super::res_grp_builder::ResGrpBuilder;

pub struct ZoomBuilder<'a> {
    from_grp: &'a ResGroup,
    new_res: u32,
    n_new_bins: usize,
    chunksize: usize,
    bin_table: (Array1<u32>, Array1<u64>, Array1<u64>),
    tig_offsets: Array1<u32>
}

impl<'a> ZoomBuilder<'a> {
    pub fn new(from_grp: &'a ResGroup, tig_lengths: ArrayView1<u64>, new_res: u32, chunksize: usize) -> ZoomBuilder<'a> {
        let tig_offsets = ZoomBuilder::build_tig_offsets(new_res, tig_lengths);
        let n_new_bins = if !tig_offsets.is_empty() {tig_offsets[tig_offsets.len() - 1] as usize} else {0};
        let bin_table = ZoomBuilder::build_bin_table_from_lengths(n_new_bins, new_res as u64, tig_lengths);

        ZoomBuilder {
            from_grp,
            new_res,
            n_new_bins,
            chunksize,
            bin_table,
            tig_offsets,
        }
    }
}

impl<'a> ResGrpBuilder for ZoomBuilder<'a> {
    fn get_resolution(&self) -> u32 {
        self.new_res
    }

    fn get_tig_offsets_view(&self) -> ArrayView1<u32> {
        self.tig_offsets.view()
    }

    fn get_bin_table(&self) -> (ArrayView1<u32>, ArrayView1<u64>, ArrayView1<u64>) {
        (self.bin_table.0.view(), self.bin_table.1.view(), self.bin_table.2.view())
    }

    fn get_bin_offsets(&self, pixels: &[PixelT]) -> Array1<u32> {
        ZoomBuilder::build_bin_offsets_from_pixels(self.n_new_bins, pixels)
    }

    fn get_pixels(&self) -> Result<Vec<PixelT>, Box<dyn Error>> {
        let new_res = self.new_res as u32;
        let bscs: Array1<(u32, u32)> = self.from_grp.get_bin_coords()?;
        let mut pixels:AHashMap<(u32, u32), u32> = AHashMap::default();

        for chunk in self.from_grp.get_raw_pixel_iter(self.chunksize) {
            let (bins1, bins2, counts) = chunk;
            for (i, &count) in counts.iter().enumerate() {
                let bin1 = bins1[i] as usize;
                let bin2 = bins2[i] as usize;
                let (crom_id1, anchor1) = bscs[bin1];
                let (crom_id2, anchor2) = bscs[bin2];
                let offset1 = self.tig_offsets[crom_id1 as usize];
                let offset2 = self.tig_offsets[crom_id2 as usize];
                let new_bin1_id = offset1 + (anchor1 / new_res) as u32;
                let new_bin2_id = offset2 + (anchor2 / new_res) as u32;
                assert!(new_bin1_id <= new_bin2_id);

                let c = pixels.entry((new_bin1_id, new_bin2_id)).or_insert(0);
                *c += count;
            }
        }

        let mut pixels: Vec<PixelT> = pixels.into_iter().map(|x| ((x.0).0, (x.0).1, x.1)).collect();
        pixels.sort_by_key(|rec| { (rec.0, rec.1) });
        Ok(pixels)
    }
}

// for (bin1, bin2, count) in self.from_grp.get_pixels()?.view() {
// let (bins1, bins2, counts) = self.from_grp.get_raw_pixels()?;
// for (i, &count) in counts.iter().enumerate() {
//     let bin1 = bins1[i] as usize;
//     let bin2 = bins2[i] as usize;
//     let (crom_id1, anchor1) = bscs[bin1];
//     let (crom_id2, anchor2) = bscs[bin2];
//     let offset1 = self.tig_offsets[crom_id1 as usize];
//     let offset2 = self.tig_offsets[crom_id2 as usize];
//     let new_bin1_id = offset1 + (anchor1 / new_res) as u32;
//     let new_bin2_id = offset2 + (anchor2 / new_res) as u32;
//     assert!(new_bin1_id <= new_bin2_id);
//
//     let c = pixels.entry((new_bin1_id, new_bin2_id)).or_insert(0);
//     *c += count;
// }
