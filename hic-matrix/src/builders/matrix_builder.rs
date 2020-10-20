use std::path::Path;
use ascii::AsciiString;
use std::error::Error;

use super::super::{Matrix, utils};
use super::pair_builder::PairsBuilder;
use super::super::writer::MatrixWriter;

pub struct MatrixBuilder {}

impl MatrixBuilder {
    pub fn build_from_pairs(pairs_file: &Path, matrix_file: &Path,
                            ord_tig_lengths: &[(AsciiString, u64)],
                            rsltn: u32,
                            is_balance: bool
    ) -> Result<Matrix, Box<dyn Error>> {
        {
            let writer = MatrixWriter::new_in_writing_mode(matrix_file)?;
            let builder = PairsBuilder::new(pairs_file, ord_tig_lengths, rsltn);
            writer.write_matrix(&builder)?;
        }
        let matrix = Matrix::from_hdf_file(matrix_file)?;
        if is_balance { matrix.balance(rsltn)?; }
        Ok(matrix)
    }

    pub fn build_from_pairs_multi_res(pairs_file: &Path, matrix_file: &Path,
                                      ord_tig_lengths: &[(AsciiString, u64)],
                                      new_ress: &[u32],
                                      is_balance: bool
    ) -> Result<Matrix, Box<dyn Error>> {
        let init_matrix = MatrixBuilder::build_from_pairs(pairs_file, matrix_file, &ord_tig_lengths,
                                                            new_ress[0], is_balance)?;
        let final_matrix = MatrixBuilder::zooming(init_matrix, new_ress, is_balance)?;
        Ok(final_matrix)
    }

    fn zooming(mut matrix: Matrix, new_ress: &[u32], is_balance: bool)
        -> Result<Matrix, Box<dyn Error>> {
        let pred = utils::get_zooming_order(new_ress);

        for (&new_res, &pred_ind) in new_ress.iter().zip(pred.iter()) {
            if pred_ind != -1 {
                let prev_res = new_ress[pred_ind as usize];
                matrix.zoom(prev_res, new_res)?;
                if is_balance { matrix.balance(new_res)?; }
            }
        }

        Ok(matrix)
    }
}