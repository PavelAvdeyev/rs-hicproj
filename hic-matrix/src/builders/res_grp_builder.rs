use ndarray::{Array1, ArrayView1};
use itertools::Itertools;
use std::error::Error;

use super::super::reader::PixelT;


pub trait ResGrpBuilder {
    fn get_resolution(&self) -> u32;

    fn get_tig_offsets_view(&self) -> ArrayView1<u32>;

    fn get_bin_table(&self) -> (ArrayView1<u32>, ArrayView1<u64>, ArrayView1<u64>);

    fn get_bin_offsets(&self, pixels: &[PixelT]) -> Array1<u32>;

    fn get_pixels(&self) -> Result<Vec<PixelT>, Box<dyn Error>>;

    fn build_tig_offsets(rsltn: u32, tig_lengths: ArrayView1<u64>) -> Array1<u32> {
        let mut count = 0_u32;
        let mut tig_offsets: Array1<u32> = Array1::default(tig_lengths.len() + 1);
        for (i, len) in tig_lengths.iter().enumerate() {
            tig_offsets[i] = count;
            let n_bins = ((*len as f64) / (rsltn as f64)).ceil() as u32;
            count += n_bins;
        }
        tig_offsets[tig_lengths.len()] = count;
        tig_offsets
    }

    fn build_bin_table_from_lengths(n_bins: usize, rsltn: u64, tig_lengths: ArrayView1<u64>)
                                    -> (Array1<u32>, Array1<u64>, Array1<u64>) {
        let mut ind = 0_usize;
        let mut chrs = Array1::<u32>::default(n_bins);
        let mut starts = Array1::<u64>::default(n_bins);
        let mut ends = Array1::<u64>::default(n_bins);

        for (i, length) in tig_lengths.iter().enumerate() {
            let n_bins = ((*length as f64) / (rsltn as f64)).ceil() as u64;

            for (prev, next) in (0..=n_bins)
                .map(|x| { if x != n_bins { x * rsltn } else { *length } })
                .tuple_windows() {
                chrs[ind] = i as u32;
                starts[ind] = prev;
                ends[ind] = next;
                ind += 1;
            }
        }

        (chrs, starts, ends)
    }

    fn build_bin_offsets_from_pixels(total_bins: usize, pixels: &[PixelT]) -> Array1<u32> {
        let tb = total_bins + 1;
        let mut bin_offsets = Array1::<u32>::default(tb);
        let start_ind;

        match pixels.first() {
            Some(bin) => {
                let ind = bin.0 as usize;
                bin_offsets[ind] = 0;
                start_ind = ind;
            },
            None => return bin_offsets,
        };

        pixels.iter().enumerate().tuple_windows().for_each(|((_, prev_bin), (i, cur_bin))| {
            if prev_bin.0 != cur_bin.0 { bin_offsets[cur_bin.0 as usize] = i as u32; }
        });
        bin_offsets[tb - 1] = pixels.len() as u32;

        for i in ((start_ind + 1)..bin_offsets.len()).rev() {
            if bin_offsets[i] == 0 { bin_offsets[i] = bin_offsets[i + 1] };
        }

        bin_offsets
    }
}

// fn calc_n_bins(rsltn: u32, tig_lengths: ArrayView1<u64>) -> usize {
//     tig_lengths.iter().map(|len| { ((*len as f64) / (rsltn as f64)).ceil() as usize }).sum()
// }
