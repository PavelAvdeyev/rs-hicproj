use std::path::Path;
use ascii::AsciiString;
use std::error::Error;

use super::super::{Matrix, utils};
use super::pair_builder::PairsBuilder;
use super::super::writer::MatrixWriter;
use super::super::balancer::Strategy;

pub fn build_from_pairs(pairs_file: &Path, matrix_file: &Path,
                        ord_tig_lengths: &[(AsciiString, u64)],
                        resolution: u32
) -> Result<(), Box<dyn Error>> {
    let writer = MatrixWriter::new_in_writing_mode(matrix_file)?;
    let builder = PairsBuilder::new(pairs_file, ord_tig_lengths, resolution);
    writer.write_matrix(&builder)?;
    Ok(())
}

pub fn build_from_pairs_with_balancing(pairs_file: &Path, matrix_file: &Path,
                                       ord_tig_lengths: &[(AsciiString, u64)],
                                       resolution: u32,
                                       balance_strategy: &Strategy
) -> Result<(), Box<dyn Error>> {
    build_from_pairs(pairs_file, matrix_file, ord_tig_lengths, resolution)?;
    balance(matrix_file, resolution, balance_strategy)?;
    Ok(())
}

pub fn build_from_pairs_multi_res(pairs_file: &Path, matrix_file: &Path,
                                  ord_tig_lengths: &[(AsciiString, u64)],
                                  resolutions: &[u32],
                                  balance_strategy: &Strategy
) -> Result<(), Box<dyn Error>> {
    build_from_pairs_with_balancing(pairs_file, matrix_file, &ord_tig_lengths,
                                                   resolutions[0], balance_strategy)?;
    zoom_with_balancing(matrix_file, &resolutions[1..], balance_strategy)?;
    Ok(())
}

pub fn balance(matrix_file: &Path, rstln: u32, balance_strategy: &Strategy) -> Result<Matrix, Box<dyn Error>> {
    let matrix = Matrix::from_hdf_file(matrix_file)?;
    matrix.balance(rstln, balance_strategy)?;
    Ok(matrix)
}

pub fn zoom_with_balancing(matrix_file: &Path, new_resolutions: &[u32],
                           balance_strategy: &Strategy) -> Result<Matrix, Box<dyn Error>> {
    let mut matrix = Matrix::from_hdf_file(matrix_file)?;

    let mut rsltns = matrix.get_resolutions();  rsltns.extend(new_resolutions);
    let pred = utils::get_zooming_order(&rsltns);

    for (&new_res, &pred_ind) in rsltns.iter().zip(pred.iter()) {
        if pred_ind != -1 {
            let prev_res = rsltns[pred_ind as usize];
            matrix.zoom(prev_res, new_res)?;
        }
    }

    for &rstln in rsltns.iter() { matrix.balance(rstln, balance_strategy)?; }

    Ok(matrix)
}
