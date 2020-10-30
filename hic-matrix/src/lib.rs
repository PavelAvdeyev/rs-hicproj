pub mod res_group;
pub mod matrix;
pub mod writer;
pub mod reader;
pub mod errors;

mod builders;
mod selector;
mod utils;
mod balancer;

use std::path::Path;
use std::error::Error;
use self::builders::matrix_builder;


pub use self::res_group::ResGroup;
pub use self::matrix::Matrix;
pub use self::balancer::Strategy;



// pub fn create_matrix_from_pairs(pairs_file: &Path, matrix_file: &Path, ordered_tig_lengths: &[(AsciiString, u64)], resolution: u32)
//     -> Result<ResGroup, Box<dyn Error>> {
//     {
//         let builder = PairsBuilder::from_pairs(matrix_file, pairs_file, resolution, ordered_tig_lengths);
//         let writer = ResGrpWriter::from_builder(&builder)?;
//         writer.write()?;
//     }
//     let matrix = ResGroup::from_hdf_file(resolution, matrix_file)?;
//     balance(&matrix)?;
//     Ok(matrix)
// }

pub fn create_matrix_from_pairs(pairs_file: &Path, tig_length_file: &Path,
                                matrix_file: &Path, rslns: &[u32],
                                strategy: &Strategy) -> Result<(), Box<dyn Error>> {
    let ord_tig_lengths = utils::parse_tig_lengths(tig_length_file)?;
    matrix_builder::build_from_pairs_multi_res(pairs_file, matrix_file, &ord_tig_lengths, rslns, strategy)?;
    Ok(())
}


pub use self::builders::matrix_builder::balance;

pub use self::builders::matrix_builder::zoom;




// pub fn create_matrix_by_zooming<'a>(prev_mat: &'a ResGroup, ordered_tig_lengths: &[(AsciiString, u64)], new_res: u32) -> Result<ResGroup, Box<dyn Error>> {
//     // if new_res % prev_mat.get_bin_size() != 0 {
//     //     println!("New resolution is not dividable by prev bin size."); // TODO
//     //     return Ok(());
//     // }
//     println!("Zooming hic_matrix from {} to {}", prev_mat.get_bin_size(), new_res);
//     let builder = ZoomBuilder::new(prev_mat, ordered_tig_lengths, new_res);
//     let writer = ResGrpWriter::from_zoomifier(&builder)?;
//     writer.write()?;
//     let matrix = ResGroup::from_hdf_file(new_res, prev_mat.get_matrix_path())?;
//     balance(&matrix)?;
//     Ok(matrix)
// }
//
// pub fn create_multi_matrix_by_zooming(prev_mat: ResGroup, ordered_tig_lengths: &[(AsciiString, u64)], new_ress: &[u32])
//                                           -> Result<ResGroup, Box<dyn Error>> {
//     let matrix_file = prev_mat.get_matrix_path();
//     let mut resolutions = Vec::new();
//     resolutions.push(prev_mat.get_bin_size());
//     resolutions.extend_from_slice(&new_ress);
//     let pred = utils::get_order_of_zooming(&resolutions);
//
//     for (&new_res, pred_ind) in resolutions.iter().zip(pred.iter()) {
//         if *pred_ind != -1 {
//             let pred_ind = *pred_ind as usize;
//             let prev_bin_size = resolutions[pred_ind];
//             let prev_mat = ResGroup::from_hdf_file(prev_bin_size, matrix_file)?;
//             create_matrix_by_zooming(&prev_mat, ordered_tig_lengths, new_res)?;
//         }
//     }
//
//     Ok(prev_mat)
// }
