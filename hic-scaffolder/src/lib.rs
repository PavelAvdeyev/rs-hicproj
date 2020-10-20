pub mod hic_graph;
pub mod trans_updater;

use std::iter::FromIterator;
use ascii::AsciiString;
use std::collections::VecDeque;
use ahash::AHashSet;
use std::error;

use hic_matrix;
use gfa_graph::graph;

use hic_graph::HiCGraphEnsemble;

pub const TIG_LEN_CUTOFF: u64 = 100_000;
pub const PRECISION: f64 = 0.0000001;
pub const MAXFINDER_CHUNKSIZE: usize = 30_000_000;


pub fn update_matrix_with_max_trans_vals(matrix: &hic_matrix::Matrix) -> Result<(), Box<dyn error::Error>> {
    trans_updater::write_best_trans_weights(matrix, TIG_LEN_CUTOFF, MAXFINDER_CHUNKSIZE)
}

pub struct PathFinder<'a> {
    ovp_graph: graph::GFAGraph,
    hic_graph: HiCGraphEnsemble<'a>,
    paths: Vec<Vec<String>>
}

impl<'a> PathFinder<'a> {
    // impl PathFinder {
    pub fn new(graph: graph::GFAGraph, matrix: &'a hic_matrix::Matrix) -> Result<PathFinder<'a>, Box<dyn error::Error>> {
        // pub fn new(graph: graph::GFAGraph, matrix: &hic_matrix::MultiMatrix) -> hdf5::Result<PathFinder> {
        Ok(PathFinder {
            ovp_graph: graph,
            hic_graph: HiCGraphEnsemble::new(matrix, TIG_LEN_CUTOFF)?,
            paths: Vec::new()
        })
    }

    fn sort_vertices_by_length(&self) -> Vec<(AsciiString, u64)> {
        let mut vertices = Vec::from_iter(self.ovp_graph.node_names()
            .map(|nn| {
                (nn.clone(), self.ovp_graph.get_tig_length(nn.as_str()).unwrap())
            }));
        vertices.sort_by_key(|(_, l)| *l);
        vertices
    }

    fn find_preferable_orientation_wrt_graph(&self) {

    }

    fn find_next_vertex(&self, cur_v: &AsciiString, psv: &mut AHashSet<AsciiString>)
                        -> Option<AsciiString> {
        println!("Trying to find next vertex for tig");
        None
    }

    pub fn find_paths(&self) {
        let vertices = self.sort_vertices_by_length();

        let mut queue = VecDeque::new();
        let (next_vertex, _) = vertices.last().unwrap();
        queue.push_back(next_vertex.clone());

        let cur_v: AsciiString;
        // let previously_suggested_vertices = AHashSet::new();
        while !queue.is_empty() {
            let cur_v = queue.pop_front().unwrap();
            println!("Starting work with vertex {}", cur_v);
            // let sug_vs = self.hic_graph.find_best_weighted_neighbors(&cur_v, true);
            //
            // println!("{}", sug_vs.len());
            // for x in sug_vs.iter() {
            //     println!("{} ", x);
            // }
            break;
        }

    }
}

// pub mod cigar;
