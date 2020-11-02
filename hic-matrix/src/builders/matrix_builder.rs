use std::path::Path;
use ascii::AsciiString;
use std::error::Error;
use ahash::AHashSet;

use super::super::Matrix;
use super::pair_builder::PairsBuilder;
use super::super::writer::MatrixWriter;
use super::super::balancer::Strategy;
use std::iter::FromIterator;
use itertools::Itertools;

// pub fn build_from_pairs(pairs_file: &Path, matrix_file: &Path,
//                         ord_tig_lengths: &[(AsciiString, u64)],
//                         resolution: u32
// ) -> Result<(), Box<dyn Error>> {
//     let writer = MatrixWriter::new_in_writing_mode(matrix_file)?;
//     let builder = PairsBuilder::new(pairs_file, ord_tig_lengths, resolution);
//     writer.write_matrix(&builder)?;
//     Ok(())
// }


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

    let mut rsltns = get_zooming_order(&matrix.get_resolutions(), new_rslns);
    assert!(!rsltns.is_empty());

    for res in rsltns.iter() {
        if res.1 != -1 {
            let prev_res = rsltns[res.1 as usize].0;
            matrix.zoom(prev_res, res.0)?;
        }
    }

    Ok(matrix)
}


fn get_zooming_order(resolutions: &[u32], new_resolutions: &[u32]) -> Vec<(u32, i32)> {
    if resolutions.is_empty() || new_resolutions.is_empty() { return vec![]; }

    let mut rsltns = Vec::from_iter (
        resolutions.iter().chain(new_resolutions.iter()).cloned()
    );
    rsltns.sort();

    for (&pre, &cur) in rsltns.iter().tuple_windows() { if pre == cur { return vec![] } }

    let mut pred = vec![-1; rsltns.len()];
    let new_res_set: AHashSet<u32> = AHashSet::from_iter(new_resolutions.iter().cloned() );

    if new_res_set.contains(&rsltns[0]) { return vec![] }

    for (i, res) in rsltns.iter().enumerate().skip(1) {
        if new_res_set.contains(res) {
            let mut p = i as i32 - 1;

            while p >= 0 {
                if res % rsltns[p as usize] == 0 {
                    pred[i] = p;
                    break;
                }
                p -= 1;
            }
        }
    }

    rsltns.into_iter().zip(pred.into_iter()).collect()
}
