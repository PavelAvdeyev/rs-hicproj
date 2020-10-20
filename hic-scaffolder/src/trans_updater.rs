use ndarray::{Array1, ArrayView1};
use std::error::Error;
use std::path::Path;

use crate::hic_matrix::{Matrix, ResGroup};
use crate::hic_matrix::writer::{MatrixWriter, self};
use crate::hic_matrix::reader::{MatrixReader, self};


pub fn read_best_trans_weights(file_path: &Path, rstln: u32) -> hdf5::Result<Array1<f64>> {
    let reader = MatrixReader::new(file_path)?;
    let reader= reader.get_res_group_reader(rstln)?;
    let grp = reader.get_root().group("bins")?;
    let mut best_weights: Array1<f64> = reader::read_dataset(&grp, "max_val")?;
    for x in best_weights.iter_mut() {*x = if x.is_nan() { 0.0 } else { *x }}
    Ok(best_weights)
}

pub fn write_best_trans_weights(matrix: &Matrix, length_cutoff: u64, chunksize: usize) -> Result<(), Box<dyn Error>> {
    let finder = MaxInRowFinder::new(length_cutoff, chunksize);
    finder.write_trans_max_in_rows(matrix)?;
    Ok(())
}


struct MaxInRowFinder {
    length_cutoff: u64,
    chunksize: usize
}

impl MaxInRowFinder {
    pub fn new(length_cutoff: u64, chunksize: usize) -> MaxInRowFinder {
        MaxInRowFinder {
            length_cutoff,
            chunksize
        }
    }

    pub fn write_trans_max_in_rows(&self, matrix: &Matrix) -> Result<(), Box<dyn Error>> {
        let tig_lengths = matrix.lengths_view();

        for rstln in matrix.get_resolutions() {
            println!("Adding max trans interaction value for each row. Resolution {}", rstln);
            let max_vals = self.calc_trans_max_in_rows(matrix.get_local_matrix(rstln).unwrap(), tig_lengths)?;
            let writer = MatrixWriter::new_in_appending_mode(matrix.get_filepath())?;
            let root = writer.get_file_handler().group(format!("resolutions/{}", rstln).as_ref())?;
            MaxInRowFinder::write_max_values_for_rows(&root, max_vals.view())?;
        }

        Ok(())
    }

    fn calc_trans_max_in_rows(&self, res_group: &ResGroup, tig_lengths: ArrayView1<u64>) -> hdf5::Result<Array1<f64>> {
        let n_bins = res_group.get_n_bins();
        let mut maxs = Array1::zeros((n_bins,));
        let mut max_updater = |id: usize, weight: f64| { if maxs[id] < weight { maxs[id] = weight; } };

        let bin_chrs = res_group.get_bin_chr_ids()?;

        for chunk in res_group.get_balanced_pixel_iter(self.chunksize)? {
            let (bins1, bins2, weights) = chunk;
            for (i, &weight) in weights.iter().enumerate() {
                let bin1 = bins1[i] as usize;
                let bin2 = bins2[i] as usize;
                let weight = if weight.is_finite() {weight} else {0.0};

                if tig_lengths[bin_chrs[bin1] as usize] >= self.length_cutoff
                    && tig_lengths[bin_chrs[bin2] as usize] >= self.length_cutoff
                    && bin_chrs[bin1] != bin_chrs[bin2] {
                    max_updater(bin1, weight);
                    max_updater(bin2, weight);
                }
            }
        }
        maxs.iter_mut().for_each(|val| if *val == 0.0 {*val = f64::NAN});
        Ok(maxs)
    }


    fn write_max_values_for_rows<T: hdf5::H5Type>(grp: &hdf5::Group, max_vals: ArrayView1<T>) -> hdf5::Result<()> {
        let grp = grp.group("bins")?;
        writer::write_dataset(&grp, "max_val", max_vals.len(), max_vals)?;
        Ok(())
    }

}

