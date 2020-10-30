use std::path::Path;
use ascii::AsciiString;
use std::error::Error;

use super::super::{Matrix, utils};
use super::pair_builder::PairsBuilder;
use super::super::writer::MatrixWriter;
use super::super::balancer::Strategy;

// pub fn build_from_pairs(pairs_file: &Path, matrix_file: &Path,
//                         ord_tig_lengths: &[(AsciiString, u64)],
//                         resolution: u32
// ) -> Result<(), Box<dyn Error>> {
//     let writer = MatrixWriter::new_in_writing_mode(matrix_file)?;
//     let builder = PairsBuilder::new(pairs_file, ord_tig_lengths, resolution);
//     writer.write_matrix(&builder)?;
//     Ok(())
// }

//TODO rewrite zoom determination

pub fn build_from_pairs(pairs_file: &Path, matrix_file: &Path,
                                       ord_tig_lengths: &[(AsciiString, u64)],
                                       resolution: u32,
                                       strategy: &Strategy
) -> Result<(), Box<dyn Error>> {
    let writer = MatrixWriter::new_in_writing_mode(matrix_file)?;
    let builder = PairsBuilder::new(pairs_file, ord_tig_lengths, resolution);
    writer.write_matrix(&builder)?;
    balance(matrix_file, &vec![resolution], strategy)?;
    Ok(())
}

pub fn build_from_pairs_multi_res(pairs_file: &Path, matrix_file: &Path,
                                  ord_tig_lengths: &[(AsciiString, u64)],
                                  rslns: &[u32],
                                  strategy: &Strategy
) -> Result<(), Box<dyn Error>> {
    build_from_pairs(pairs_file, matrix_file, &ord_tig_lengths, rslns[0], strategy)?;
    zoom(matrix_file, &rslns[1..])?;
    balance(matrix_file, &rslns[1..], strategy)?;
    Ok(())
}

pub fn balance(matrix_file: &Path, rslns: &[u32], strategy: &Strategy) -> Result<Matrix, Box<dyn Error>> {
    let matrix = Matrix::from_hdf_file(matrix_file)?;
    for &r in rslns { matrix.balance(r, strategy)? };
    Ok(matrix)
}

pub fn zoom(matrix_file: &Path, new_rslns: &[u32]) -> Result<Matrix, Box<dyn Error>> {
    let mut matrix = Matrix::from_hdf_file(matrix_file)?;

    let mut rsltns = matrix.get_resolutions();  rsltns.extend(new_rslns);
    let pred = utils::get_zooming_order(&rsltns);

    for (&new_res, &pred_ind) in rsltns.iter().zip(pred.iter()) {
        if pred_ind != -1 {
            let prev_res = rsltns[pred_ind as usize];
            matrix.zoom(prev_res, new_res)?;
        }
    }

    Ok(matrix)
}
