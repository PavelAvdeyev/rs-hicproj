use std::path::{Path, PathBuf};
use ahash::AHashMap;
use ndarray::{Array1, ArrayView1};
use std::iter::FromIterator;
use ascii::{AsciiString, AsciiStr};
use std::error::Error;

use super::res_group::ResGroup;
use super::reader::MatrixReader;
use super::balancer::{Balancer, Strategy};
use super::writer::MatrixWriter;
use super::builders::zoom_builder::ZoomBuilder;
use super::errors::MatrixResolutionError;

const ZOOM_CHUNKSIZE: usize = 30_000_000;

#[derive(Default,Debug)]
pub struct Matrix {
    resolutions: AHashMap<u32, ResGroup>,
    name2order: AHashMap<AsciiString, usize>,
    tig_order: Array1<AsciiString>,
    tig_lengths: Array1<u64>,
    file_path: PathBuf,
}


impl<'a> Matrix {
    pub fn new() -> Matrix {
        Matrix {
            resolutions: AHashMap::default(),
            name2order: Default::default(),
            tig_order: Default::default(),
            tig_lengths: Default::default(),
            file_path: Default::default()
        }
    }

    pub fn from_hdf_file(file_path: &Path) -> hdf5::Result<Matrix> {
        let reader = MatrixReader::new(file_path)?;
        let (tig_order, tig_lengths) = reader.read_chroms_info()?;
        let mut matrix = Matrix {
            resolutions: AHashMap::default(),
            name2order: tig_order.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect(),
            tig_order,
            tig_lengths,
            file_path: PathBuf::from(file_path)
        };

        let resolutions = reader.read_resolutions()?;
        for res in resolutions.into_iter() {
            matrix.register_new_resolution(res)?;
        }

        Ok(matrix)
    }

    pub fn init_selectors(mut self) -> hdf5::Result<Matrix> {
        for (_, m) in self.resolutions.iter_mut() {
            m.init_selector()?;
        }
        Ok(self)
    }

    pub fn balance(&self, rstln: u32, strategy: &Strategy) -> Result<(), Box<dyn Error>> {
        println!("Balance {}", rstln);
        match self.resolutions.get(&rstln) {
            Some(res_group) => {
                let balancer = Balancer::new();
                let weights = match strategy {
                    Strategy::ICGenomeWide => balancer.balance_by_ic_genomewide(res_group),
                    Strategy::BinLength => Some(balancer.balance_by_resolution(res_group)),
                    Strategy::None => None
                };

                if let Some(wghs) = weights {
                    let writer = MatrixWriter::new_in_appending_mode(self.file_path.as_path())?;
                    writer.write_balancing_weights(rstln, wghs.view())?;
                }

                Ok(())
            },
            _ => Err(MatrixResolutionError.into())
        }
    }

    pub fn zoom(&mut self, from_rstln: u32, to_rstln: u32) -> Result<(), Box<dyn Error>> {
        println!("Zooming matrix from {} to {}", from_rstln, to_rstln);
        match self.resolutions.get(&from_rstln) {
            Some(from_grp) => {
                {
                    let builder = ZoomBuilder::new(from_grp, self.tig_lengths.view(), to_rstln, ZOOM_CHUNKSIZE);
                    let writer = MatrixWriter::new_in_appending_mode(self.file_path.as_path())?;
                    writer.write_resolution_group(&builder)?;
                }
                self.register_new_resolution(to_rstln)?;
                Ok(())
            },
            _ => Err(MatrixResolutionError.into())
        }
    }

    pub fn get_filepath(&self) -> &Path {
        self.file_path.as_path()
    }

    pub fn get_n_chroms(&self) -> usize {
        self.tig_order.len()
    }

    pub fn get_tig_id(&self, tig: &AsciiStr) -> Option<usize> {
        self.name2order.get(tig).cloned()
    }

    pub fn get_tig_name(&self, tig_id: usize) -> Option<AsciiString> {
        if tig_id >= self.tig_order.len() { None } else { Some(self.tig_order[tig_id].clone()) }
    }

    pub fn get_local_matrix(&self, resolution: u32) -> Option<&ResGroup> {
        self.resolutions.get(&resolution)
    }

    pub fn tig_order_view(&self) -> ArrayView1<AsciiString> {
        self.tig_order.view()
    }

    pub fn lengths_view(&self) -> ArrayView1<u64> {
        self.tig_lengths.view()
    }

    pub fn get_resolutions(&self) -> Vec<u32> {
        Vec::from_iter(self.resolutions.keys().copied())
    }

    fn register_new_resolution(&mut self, rstln: u32) -> hdf5::Result<()> {
        println!("We are registering new resolution {}", rstln);
        let reader = MatrixReader::new(self.file_path.as_path())?;
        let res_group_reader = reader.get_res_group_reader(rstln)?;
        let res_group = ResGroup::new(rstln, res_group_reader)?;
        self.resolutions.insert(rstln, res_group);
        Ok(())
    }
}

//    pub fn balance_all(&self) -> Result<(), Box<dyn Error>> {
//         println!("Balance all resolutions!");
//         for &res in self.resolutions.keys() {
//             self.balance(res)?;
//         }
//         Ok(())
//     }