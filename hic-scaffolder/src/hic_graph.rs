use ascii::AsciiString;
use ahash::{AHashMap, AHashSet};
use ndarray::{Array1, ArrayView1};
use std::error;

use super::trans_updater;
use super::PRECISION;
// use crate::hic_matrix;
// use crate::hic-scaffolder::PRECISION;
// use crate::hic_scaffolder::trans_updater;


pub struct HiCGraphEnsemble<'a> {
    graphs: Vec<HiCGraph<'a>>,
    resolutions: Vec<u32>,
    matrix: &'a hic_matrix::Matrix,
}

impl<'a> HiCGraphEnsemble<'a> {
    pub fn new(matrix: &'a hic_matrix::Matrix, length_cutoff: u64) -> Result<HiCGraphEnsemble<'a>, Box<dyn error::Error>> {
        let mut resolutions = matrix.get_resolutions(); resolutions.sort();
        let mut graphs = Vec::new();
        for &res in resolutions.iter() {
            graphs.push(HiCGraph::new(matrix.get_local_matrix(res).unwrap(),
             trans_updater::read_best_trans_weights(matrix.get_filepath(), res)?,
               matrix.lengths_view(), length_cutoff)?
            );
        }

        Ok(HiCGraphEnsemble {
            graphs,
            resolutions,
            matrix,
        })
    }

    pub fn find_best_weighted_neighbors(&self, cur_v: &AsciiString, is_start: bool) -> AHashSet<AsciiString> {
        let mut svs = AHashSet::new();
        for (i, graph) in self.graphs.iter().enumerate() {
            println!("Working with resolution {}", self.resolutions[i]);
            if let Some(tig_id) = self.matrix.get_tig_id(cur_v) {
                println!("Our tig id is {}", tig_id);
                let o_bwn = graph.find_best_weighted_neighbor(tig_id as u32, is_start);
                if let Some(bwn) = o_bwn {
                    let next_tig_name = self.matrix.get_tig_name(bwn.0 as usize).expect("Something terrible happened. ");
                    println!("Found the best neighbor {} with weight {}", next_tig_name, bwn.1);
                    svs.insert(next_tig_name);
                    // break;
                }
            }
        }
        svs
    }
}

struct HiCGraph<'a> {
    matrix: &'a hic_matrix::ResGroup,
    start_bins: AHashMap<u32, u32>,
    end_bins: AHashMap<u32, u32>,
    bin_info: Array1<(u32, u32)>, //chr_id, start_t
    tig_lengths: ArrayView1<'a, u64>,
    best_weights: Array1<f64>,
    length_cutoff: u64
}

impl<'a> HiCGraph<'a> {
    pub fn new(matrix: &'a hic_matrix::ResGroup,
               best_weights: Array1<f64>,
               tig_lengths: ArrayView1<'a, u64>,
               length_cutoff: u64)
        -> Result<HiCGraph<'a>, Box<dyn error::Error>> {
        println!("Building graph that corresponds to resolution {}", matrix.get_resolution());
        let mut graph = HiCGraph {
            matrix,
            start_bins: Default::default(),
            end_bins: Default::default(),
            bin_info: Default::default(),
            tig_lengths,
            best_weights,
            length_cutoff
        };
        graph.init_bin_info()?;
        Ok(graph)
    }

    pub fn find_best_weighted_neighbor(&self, name: u32, is_start: bool) -> Option<(u32, f64)> {
        let bin_id = if is_start { self.start_bins.get(&name) } else { self.end_bins.get(&name) };

        bin_id.and_then(|id| {
            let (row_chr, _) = self.bin_info[*id as usize];
            // println!("We are working with row {} and chr {}", id, row_chr);

            let mut max_elem = None;
            let row = self.matrix.get_balanced_row_as_nnz_elems(*id as usize)
                .expect("Problem with requesting row from matrix.");

            for (cid, cw) in row.iter() {
                let (cchr, _) = self.bin_info[*cid as usize];

                // println!("We are looking at cell {},{} and chrs {},{}", id, cid, row_chr, cchr);
                if row_chr == cchr {continue};
                if self.tig_lengths[cchr as usize] < self.length_cutoff {continue;}

                let o_bbw = self.get_best_buddy_weight(*id, *cid, *cw);

                if let Some(bbw) = o_bbw {
                    if bbw - 1.0 < PRECISION { continue }

                    max_elem = match max_elem {
                        Some((_, me)) => if bbw > me {Some((cchr, bbw))} else {max_elem},
                        None => Some((cchr, bbw))
                    }
                }
            }
            max_elem
        })
    }

    fn get_best_buddy_weight(&self, bin_id1: u32, bin_id2: u32, cur_weight: f64) -> Option<f64> {
        // println!("Trying to find best weight for {} {}", bin_id1, bin_id2);
        let max1 = self.best_weights[bin_id1 as usize];
        let max2 = self.best_weights[bin_id2 as usize];
        // println!("Best weights in corresponding rows are {} {} and current weight {}", max1, max2, cur_weight);

        assert!(max1.is_finite());
        assert!(max2.is_finite());
        let max_w = if max1 > max2 {max1} else {max2};

        if max_w == 0.0 {
            None
        } else if max_w - cur_weight < PRECISION {
            // println!("We need to find new weight");
            let max_w = self.get_max_weight(bin_id1, bin_id2);
            // println!("We found new max weight {:?}", max_w);
            if let Some(w) = max_w { Some(cur_weight / w) } else { None }
        } else {
            Some(cur_weight / max_w)
        }
    }

    fn get_max_weight(&self, bin_id1: u32, bin_id2: u32) -> Option<f64> {
        let mut chrs = AHashSet::<u32>::new();
        chrs.insert(self.bin_info[bin_id1 as usize].0);
        chrs.insert(self.bin_info[bin_id2 as usize].0);

        let row = self.matrix.get_balanced_row_as_nnz_elems(bin_id1 as usize)
                            .expect("Problem with requesting row from matrix.");
        let max1 = self.find_max_elem_in_row_wrt_chrs(&row, &chrs);

        let row = self.matrix.get_balanced_row_as_nnz_elems(bin_id2 as usize)
            .expect("Problem with requesting row from matrix.");
        let max2 = self.find_max_elem_in_row_wrt_chrs(&row, &chrs);

        max1.zip(max2).map(|(a1, b1)| {
            if a1.1 > b1.1 {a1.1} else {b1.1}
        })
    }

    fn find_max_elem_in_row_wrt_chrs(&self, row: &[(u32, f64)], chrs: &AHashSet<u32>) -> Option<(u32, f64)> {
        let mut max_elem = None;

        for (i, v) in row.iter() {
            let (c_chr, _) = self.bin_info[*i as usize];
            if self.tig_lengths[c_chr as usize] < self.length_cutoff {continue;}
            if chrs.get(&c_chr).is_none() {
                max_elem = match max_elem {
                    Some((_, v_m)) => if *v > v_m { Some((*i, *v)) } else { max_elem },
                    None => Some((*i, *v)),
                };
            }
        }

        max_elem
    }

    fn init_bin_info(&mut self) -> hdf5::Result<()> {
        self.bin_info = self.matrix.get_bin_coords()?;

        let tig_offsets = self.matrix.get_tigs_offsets()?;
        for i in 0..(tig_offsets.len() - 1) {
            let (start_bin, end_bin) = (tig_offsets[i], tig_offsets[i + 1] - 1);
            let (chr, _) = self.bin_info[start_bin as usize];
            self.start_bins.insert(chr, start_bin);
            self.end_bins.insert(chr, end_bin);
        }

        Ok(())
    }

}

//
// let mut chrs = AHashSet::<i32>::new();
// chrs.insert(self.bin_info[bin_id1 as usize].0);
// chrs.insert(self.bin_info[bin_id2 as usize].0);
//
// match max1.zip(max2) {
//     Some((a, b)) => {
//         let mut max_t = if a.1 > b.1 {a} else {b};
//
//         if max_t.1 - cur_weight < PRECISION {
//             let row = self.matrix.get_balanced_row_as_nnz_elems(bin_id1 as usize)
//                 .expect("Problem with requesting row from matrix.");
//             let max1 = self.find_max_elem_in_row_wrt_chrs(&row, &chrs);
//             let row = self.matrix.get_balanced_row_as_nnz_elems(bin_id2 as usize)
//                 .expect("Problem with requesting row from matrix.");
//             let max2 = self.find_max_elem_in_row_wrt_chrs(&row, &chrs);
//
//             if let Some((a1, b1)) = max1.zip(max2) {
//                 let max_t = if a1.1 > b1.1 {a1} else {b1};
//                 Some(cur_weight / max_t.1)
//             } else {
//                 None
//             }
//         } else {
//             Some(cur_weight / max_t.1)
//         }
//     },
//     None => None,
// }