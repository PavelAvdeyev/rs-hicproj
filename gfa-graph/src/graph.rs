use std::collections::VecDeque;
use std::str::FromStr;

use petgraph::Graph;
use ahash::{AHashMap, AHashSet};
use ascii::AsciiString;
use petgraph::graph::NodeIndex;

use super::utils::Orientation;
use super::parser::prepack1::Gfa1Prepack;
use super::parser::structs1::{SegRec, LinkRec};

pub struct GFAGraph {
    grh: Graph<(), ()>,
    sequences: AHashMap<AsciiString, AsciiString>,
    seq_lengths: AHashMap<AsciiString, u64>,
    name2index: AHashMap<AsciiString, NodeIndex>, // valid until deletion
    index2name: AHashMap<NodeIndex, AsciiString>, // valid until deletion
}

impl GFAGraph {
    pub fn new() -> GFAGraph {
        GFAGraph {
            grh: Graph::<(), ()>::new(),
            sequences: AHashMap::default(),
            seq_lengths: AHashMap::default(),
            name2index: Default::default(),
            index2name: Default::default()
        }
    }

    pub fn from_prepack(prepack: &Gfa1Prepack) -> GFAGraph {
        println!("Building graph from prepack.");
        let mut graph = GFAGraph::new();

        println!("Adding nodes.");
        for rec in prepack.seq_recs_iter() {
            graph.add_node(rec);
        }

        println!("Adding links.");
        for rec in prepack.link_recs_iter() {
            graph.add_edge_from_link(rec);
        }

        graph
    }

    pub fn node_names(&self) -> impl Iterator<Item = &AsciiString> {
        self.sequences.keys()
    }

    pub fn get_tig_length(&self, name: &str) -> Option<u64> {
        AsciiString::from_str(name).ok()
            .and_then(|key| self.seq_lengths.get(&key) )
            .map(|&v| v)
    }

    // pub fn to_prepack(&self) {
    //     println!("Converting assembly graph to GFA records");
    //     let mut seq_recs: Vec<SeqRec> = Vec::from_iter(self.sequences.iter().map(|(name, seq)| {
    //
    //     }));
    // }

    pub fn has_path(&self, source: &AsciiString, target: &AsciiString) -> bool {
        let mut visited = AHashSet::default();
        let s;
        let t;

        match self.name2index.get(source).zip(self.name2index.get(target)) {
            Some((si, ti)) => {
                s = *si;
                t = *ti;
            },
            _ => return false,
        };

        fn dfs(grh: &Graph::<(), ()>, cur: NodeIndex, tgt: NodeIndex, visited: &mut AHashSet<NodeIndex>) -> bool {
            visited.insert(cur.clone());
            let mut answer = false;
            for node in grh.neighbors(cur) {
                if visited.get(&node).is_none() {
                    if node == tgt {
                        return true;
                    };
                    answer |= dfs(grh, node, tgt, visited);
                }
            }
            answer
        }
        dfs(&self.grh, s, t, &mut visited)
    }

    // def short_paths_via_bfs(self, source):
    //         q = Queue()
    //         visited = set()
    //         dists = defaultdict()
    //
    //         q.put(source)
    //         dists[source] = 0
    //
    //         while not q.empty():
    //             current_node = q.get()
    //             # logger.debug(f"Working with {current_node}")
    //             # visited.add(self.get_tig_name(current_node))
    //
    //             if self.graph.out_degree(current_node) == 0:
    //                 comp_current_node = self.get_complement_node(current_node)
    //                 if comp_current_node not in dists:
    //                     dists[comp_current_node] = dists[current_node] + 1
    //                 # logger.debug(f"Replacing {current_node} with {comp_current_node}")
    //                 current_node = comp_current_node
    //
    //             # logger.debug(f"Successors {list(self.graph.successors(current_node))}")
    //             for next_node in self.graph.successors(current_node):
    //                 # self.get_tig_name(next_node) not in visited and
    //                 if next_node not in dists:
    //                     q.put(next_node)
    //                     dists[next_node] = dists[current_node] + 1
    //
    //         return dists

    fn out_degree(&self, source: NodeIndex) -> u32 {
        self.grh.neighbors(source).fold(0_u32, |x, _| x + 1)
    }

    pub fn short_paths_via_bfs(&self, source: &AsciiString) -> AHashMap<AsciiString, i32> {
        let mut queue = VecDeque::new();
        let mut dists = AHashMap::default();

        queue.push_back(source.clone());
        dists.insert(source.clone(), 0);

        while !queue.is_empty() {
            let cur = queue.pop_back().unwrap();

            // if self.out_degree(*self.name2index.get(&cur).unwrap()) == 0 {
            //     let com_cur;
            //     cur = com_cur;
            // }

            for nn in self.grh.neighbors(*self.name2index.get(&cur).unwrap()) {
                if dists.get(&cur).is_none() {
                    let cur_dist = *dists.get(&cur).unwrap();
                    let str_name = self.index2name.get(&nn).unwrap();
                    queue.push_back(str_name.clone());
                    dists.insert(str_name.clone(),cur_dist + 1);
                }
            }
        }

        dists
    }

    fn add_edge_from_link(&mut self, rec: &LinkRec) {
        let from_node_fow = GFAGraph::get_fow_node_name(&rec.from_name);
        let from_node_rev = GFAGraph::get_rev_node_name(&rec.from_name);
        let to_node_fow = GFAGraph::get_fow_node_name(&rec.to_name);
        let to_node_rev = GFAGraph::get_rev_node_name(&rec.to_name);

        match (&rec.from_strand, &rec.to_strand) {
            (Orientation::Forward, Orientation::Forward) => {
                self.add_edge(&from_node_fow, &to_node_fow);
                self.add_edge(&to_node_rev, &from_node_rev);
            },
            (Orientation::Forward, Orientation::Reverse) => {
                self.add_edge(&from_node_fow, &to_node_rev);
                self.add_edge(&to_node_fow, &from_node_rev);
            },
            (Orientation::Reverse, Orientation::Forward) => {
                self.add_edge(&from_node_rev, &to_node_fow);
                self.add_edge(&to_node_rev, &from_node_fow);
            },
            (Orientation::Reverse, Orientation::Reverse) => {
                self.add_edge(&from_node_rev, &to_node_rev);
                self.add_edge(&to_node_fow, &from_node_fow);
            }
        }
    }

    fn add_edge(&mut self, from_node: &AsciiString, to_node: &AsciiString) {
        let from_ni = self.name2index.get(from_node);
        let to_ni = self.name2index.get(to_node);
        match from_ni.zip(to_ni) {
            Some((fni, tni)) => {
                self.grh.add_edge(*fni, *tni, ());
            }
            _ => ()
        }
    }

    fn add_node(&mut self, rec: &SegRec) {
        if let Some(seq) = rec.seq.clone() {
            self.sequences.insert(rec.name.clone(), seq);
        }

        if let Some(sl) = rec.get_length() {
            self.seq_lengths.insert(rec.name.clone(), sl);
        }

        let ni: NodeIndex = self.grh.add_node(());
        let node_name = GFAGraph::get_fow_node_name(&rec.name);
        self.name2index.insert(node_name.clone(), ni);
        self.index2name.insert(ni, node_name);

        let ni: NodeIndex = self.grh.add_node(());
        let node_name = GFAGraph::get_rev_node_name(&rec.name);
        self.name2index.insert(node_name.clone(), ni);
        self.index2name.insert(ni, node_name);
    }

    pub fn get_fow_node_name(name: &AsciiString) -> AsciiString {
        AsciiString::from(Orientation::to_node_name(name.clone(), Orientation::Forward))
    }

    pub fn get_rev_node_name(name: &AsciiString) -> AsciiString {
        AsciiString::from(Orientation::to_node_name(name.clone(), Orientation::Reverse))
    }
}