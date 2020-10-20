use std::path::{PathBuf, Path};
use std::iter::FromIterator;
use std::io::{self, BufWriter, Write};
use std::fs::File;
use log::{info, trace, warn};

use itertools::Itertools;
use ascii::AsciiString;
use bam::{RecordReader, Header};

use super::pair_record::{self, PairRecord};

// When a read matches in its entirety, with an equal score in multiple locations, one of the locations is picked at
// random, is labeled as primary, will be given a mapping quality of zero and will have an XA tag that contains the
// alternative locations (this is identical to how bwa aln worked)
//
// When different, non-overlapping regions of a read align with high scores to different, non-linear locations in
// the genome, the higher score alignment will be labeled as primary, the others may be reported as secondary alignments.
// There is some threshold on how many of these secondary alignments will be reported
//
// When two complementary regions of a read (the two pieces add up to the full read) align to two different, non-linear
// genomic locations one of the alignment will be labeled as primary, the other as supplementary alignment

const MATCHED_RATE_TRESH: f64 = 0.8;
const MIN_MAPQ: u8 = 10;
const MAX_MOLECULE_SIZE: u64 = 2000;

pub enum RescueStrategy {
    Unique,
    Simple,
    Complex
}

pub struct Converter {
    bam_path: PathBuf,
    _graph: Option<PathBuf>,
    pair_file: BufWriter<File>,
    strategy: RescueStrategy,
    max_molecule_size: u64,
    matched_rate_tresh: f64,
    min_mapq: u8,
    mapq_zero_rescue: bool,
    stats: ConverterStat
}

impl Converter {
    pub fn new(bam_file: &Path, _graph: Option<PathBuf>, pair_file: &Path) -> Converter {
        Converter {
            bam_path: PathBuf::from(bam_file),
            _graph,
            pair_file: BufWriter::new(File::create(pair_file).expect("Problem with file")),
            strategy: RescueStrategy::Complex,
            max_molecule_size: MAX_MOLECULE_SIZE,
            matched_rate_tresh: MATCHED_RATE_TRESH,
            min_mapq: MIN_MAPQ,
            mapq_zero_rescue: false,
            stats: ConverterStat::new()
        }
    }

    pub fn update_min_mapq(mut converter: Converter, mapq: u8) -> Converter {
        converter.min_mapq = mapq;
        converter
    }

    pub fn update_max_mol_size(mut converter: Converter, mol_sz: u64) -> Converter {
        converter.max_molecule_size = mol_sz;
        converter
    }

    pub fn update_matched_rate_tresh(mut converter: Converter, tresh: f64) -> Converter {
        converter.matched_rate_tresh = tresh;
        converter
    }

    pub fn update_mapq_zero_rescue(mut converter: Converter, mapq_zero_rescue: bool) -> Converter {
        converter.mapq_zero_rescue = mapq_zero_rescue;
        converter
    }

    pub fn save_statistic(&self, file_path: &Path) {
        self.stats.dump_stats_to_file(file_path);
    }

    pub fn convert(&mut self) -> io::Result<()> {
        let mut recs1 = vec![];
        let mut recs2 = vec![];
        let mut prev_read_id: Option<AsciiString> = None;

        let mut reader = bam::BamReader::from_path(self.bam_path.as_path(), 0).unwrap();

        trace!("Reading header...");
        let header = reader.header().clone();

        trace!("Reading body...");
        let mut record = bam::Record::new();
        loop {
            let is_success = reader.read_into(&mut record)?;

            if !is_success { break; }

            self.stats.update_align_count();

            let r_id = AsciiString::from_ascii(record.name()).unwrap();
            if let Some(pr_id) = &prev_read_id {
                if pr_id.as_str() != r_id.as_str() {
                    trace!("New group of alignments for {} is detected with {} and {}.", pr_id, recs1.len(), recs2.len());
                    self.parse_paired_alignments(&recs1, &recs2, &header);
                    recs1.clear(); recs2.clear();
                }
            }

            prev_read_id = Some(r_id);

            if record.flag().first_in_pair() {
                recs1.push(record.clone())
            } else {
                assert!(record.flag().last_in_pair());
                recs2.push(record.clone())
            }

            if self.stats.alignment_counter % 1000000 == 0 {
                info!("{} alignments were processed", self.stats.alignment_counter);
            }
        }

        trace!("Dump the latest group of alignments");
        self.parse_paired_alignments(&recs1, &recs2, &header);

        Ok(())
    }

    fn parse_paired_alignments(&mut self, recs1: &[bam::Record], recs2: &[bam::Record], header: &bam::Header) {
        if recs1.is_empty() || recs2.is_empty() {
            return;
        }

        self.stats.update_read_count();
        self.stats.update_alignment_count(recs1);
        self.stats.update_alignment_count(recs2);

        let prim_r1 = self.get_primary_alignment(recs1);
        let prim_r2 = self.get_primary_alignment(recs2);

        if prim_r1.zip(prim_r2).is_none() {
            warn!("It must be one and only one primary alignment for each read in pair");
            return;
        }

        let prim_r1 = prim_r1.unwrap();
        let prim_r2 = prim_r2.unwrap();

        self.stats.update_mapping_count(&prim_r1, &prim_r2);

        if !prim_r1.flag().is_mapped() || !prim_r2.flag().is_mapped() {
            trace!("At least one of reads in pair are unmapped.");
            return;
        }

        self.stats.update_mapq_count(prim_r1, prim_r2);

        if recs1.len() == 1 && recs2.len() == 1 {
            trace!("Pair read aligned 1&1 (perfectly) .");
            let hic_records = self.convert_to_pair_records(prim_r1, prim_r2, header);
            self.write_records(PairType::UU, &hic_records);
        } else if (recs1.len() == 1 || recs2.len() == 1)
            && matches!(self.strategy, RescueStrategy::Simple | RescueStrategy::Complex) {
            trace!("Pair read aligned as 1&2 (simple).");
            let resc_linear_pair = self.rescue_simple_walk(recs1, recs2);
            if let Some((rec1, rec2)) = resc_linear_pair {
                trace!("Hi-C read was rescued successfully.");
                let hic_records = self.convert_to_pair_records(rec1, rec2, header);
                self.write_records(PairType::UD, &hic_records);
            }
        } else if recs1.len() < 3 && recs2.len() < 3 && matches!(self.strategy, RescueStrategy::Complex) {
            trace!("Pair read aligned as 2&2 (complex).");
            let resc_linear_pair = self.rescue_complex_walk(recs1, recs2);
            if let Some((rec1, rec2)) = resc_linear_pair {
                trace!("Hi-C read was rescued successfully.");
                let hic_records = self.convert_to_pair_records(rec1, rec2, header);
                self.write_records(PairType::DD, &hic_records);
            }
        }
    }

    fn get_primary_alignment<'a>(&self, records: &'a [bam::Record]) -> Option<&'a bam::Record> {
        let mut primary = None;
        for rec in records {
            if !rec.flag().is_secondary() && !rec.flag().is_supplementary() {
                primary = Some(rec)
            }
        }
        primary
    }

    fn rescue_simple_walk<'a>(&self, recs1: &'a[bam::Record], recs2: &'a[bam::Record])
        -> Option<(&'a bam::Record, &'a bam::Record)> {
        if recs1.len() != 1 && recs2.len() != 1 { return None; }

        let linear_algn = if recs1.len() == 1 { &recs1[0] } else { &recs2[0] };
        let (falgn, salgn) = if recs1.len() == 1 {(&recs2[0], &recs2[1])} else {(&recs1[0], &recs1[1])};

        let dist_fa = pair_record::calc_dist(falgn, linear_algn);
        let dist_sa = pair_record::calc_dist(salgn, linear_algn);

        trace!("Distances for simple read are {} {}", dist_fa, dist_sa);
        let cor_end;
        let on_linear_side;
        if dist_sa < self.max_molecule_size {
            cor_end = falgn;
            on_linear_side = salgn;
        } else if dist_fa < self.max_molecule_size {
            cor_end = salgn;
            on_linear_side = falgn;
        } else {
            return None
        }

        if pair_record::is_opposite_pair(on_linear_side, linear_algn) {
            Some((cor_end, linear_algn))
        } else {
            None
        }
    }

    fn rescue_complex_walk<'a>(&self, recs1: &'a[bam::Record], recs2: &'a[bam::Record])
                               -> Option<(&'a bam::Record, &'a bam::Record)> {
        if recs1.len() != 2 || recs2.len() != 2 { return None; }

        let dist00 = pair_record::calc_dist(&recs1[0], &recs2[0]);
        let dist01 = pair_record::calc_dist(&recs1[0], &recs2[1]);
        let dist10 = pair_record::calc_dist(&recs1[1], &recs2[0]);
        let dist11 = pair_record::calc_dist(&recs1[1], &recs2[1]);

        trace!("Distances for complex read are {} {} {} {}", dist00, dist01, dist10, dist11);

        if dist00 < self.max_molecule_size && dist11 < self.max_molecule_size {
            if pair_record::is_opposite_pair(&recs1[0], &recs2[0])
                && pair_record::is_opposite_pair(&recs1[1], &recs2[1]) {
                Some((&recs1[0], &recs2[1]))
            } else {
                None
            }
        } else if dist01 < self.max_molecule_size && dist10 < self.max_molecule_size {
            if pair_record::is_opposite_pair(&recs1[0], &recs2[1])
                && pair_record::is_opposite_pair(&recs1[1], &recs2[0]) {
                Some((&recs1[0], &recs2[0]))
            } else {
                None
            }
        } else {
            None
        }
    }

    fn convert_to_pair_records(&self, prim_r1: &bam::Record, prim_r2: &bam::Record, header: &Header) -> Vec<PairRecord> {
        fn get_pairs(rec: &bam::Record, min_mapq: u8, is_rescue: bool) -> Vec<bam::Record> {
            let mut ans = Vec::new();

            if rec.mapq() >= min_mapq {
                ans.push(rec.clone())
            } else if rec.mapq() == 0 && is_rescue {
                // println!("Add for future support. ")
                //                 can_save = (rec.matched_proportion() - self.MATCHED_PROPORTION_TRESH >= 0)
                //
                //                 if can_save:
                //                     alts = self.graph.get_recs_within_overlap(rec.ref_name, rec.ref_algn_start)
                //
                //                     for alt in alts:
                //                         n_name, n_pos, n_strand = alt
                //                         new_rec = deepcopy(rec)
                //                         new_rec.ref_name = n_name
                //                         new_rec.ref_algn_start = n_pos
                //                         new_rec.strand = n_strand
                //                         logger.debug(f"New pair in overlaps {n_name} {n_pos} {n_strand}")
                //                         ans.append(new_rec)
                //
                //                     if len(alts):
                //                         logger.debug(f"New pair in overlaps {rec.ref_name} {rec.ref_algn_start} {rec.strand}")
                //                         ans.append(rec)
                //                         logger.debug(f"We rescued read {rec.query_name} with {len(ans)} alignments")
            }

            ans
        }

        let recs1 = get_pairs(prim_r1, self.min_mapq, self.mapq_zero_rescue);
        let recs2 = get_pairs(prim_r2, self.min_mapq, self.mapq_zero_rescue);

        Vec::from_iter(recs1.iter().cartesian_product(recs2.iter()).map(|(r1, r2)| {
            PairRecord::from_bams(r1, r2, header)
        }))
    }

    fn write_records(&mut self, tp: PairType, records: &[PairRecord]) {
        self.stats.update_cis_trans_count(records);
        self.stats.update_pair_count(tp, records.len() as u64);
        trace!("Saving {} Hi-C pairs into file", records.len());
        for rec in records {
            writeln!(self.pair_file, "{}", rec.to_string()).expect("Problem with writing file");
        }
    }
}

enum PairType {
    UU,
    UD,
    DD
}

struct ConverterStat {
    read_counter: u64,  // total number of reads
    alignment_counter: u64,  // total number of alignment records
    pairs_counter: u64, // total amount of pairst
    
    // n = non-mapped, m = mapped, each letter is for a read in pair
    nn_counter: u64, 
    nm_counter: u64,
    mm_counter: u64,

    // supp = supplementary, sec = secondary, prime = prime
    prime_counter: u64,
    supp_counter: u64,
    sec_counter: u64,

    // mapq = 0 is interesting since it may have perfect alignments but with XA tag
    mq00_counter: u64,
    mq01_counter: u64,
    mq11_counter: u64,

    // cis or trans hic pairs
    intra_counter: u64,
    inter_counter: u64,

    // pair_type of hic pairs
    uu_pair_counter: u64,
    uw_pair_counter: u64,
    ww_pair_counter: u64,
}

impl ConverterStat { 
    pub fn new() -> ConverterStat { 
        ConverterStat {
            read_counter: 0,
            alignment_counter: 0,
            pairs_counter: 0,
            nn_counter: 0,
            nm_counter: 0,
            mm_counter: 0,
            prime_counter: 0,
            supp_counter: 0,
            sec_counter: 0,
            mq00_counter: 0,
            mq01_counter: 0,
            mq11_counter: 0,
            intra_counter: 0,
            inter_counter: 0,
            uu_pair_counter: 0,
            uw_pair_counter: 0,
            ww_pair_counter: 0
        }
    }

    pub fn update_align_count(&mut self) {
        self.alignment_counter += 1;
    }

    pub fn update_read_count(&mut self) {
        self.read_counter += 1
    }

    pub fn update_pair_count(&mut self, tp: PairType, count: u64) {
        match tp {
            PairType::UU => { self.uu_pair_counter += count; },
            PairType::UD => { self.uw_pair_counter += count; }
            PairType::DD => { self.ww_pair_counter += count; }
        }
    }

    pub fn update_alignment_count(&mut self, records: &[bam::Record]) {
        for rec in records {
            if rec.flag().is_secondary() {
                self.sec_counter += 1;
            } else if rec.flag().is_supplementary() {
                self.supp_counter += 1;
            } else {
                self.prime_counter += 1;
            }
        }
    }

    pub fn update_mapping_count(&mut self, rec1: &bam::Record, rec2: &bam::Record) {
        if !rec1.flag().is_mapped() && !rec2.flag().is_mapped() {
            self.nn_counter += 1;
        } else if rec1.flag().is_mapped() != rec2.flag().is_mapped() {
            self.nm_counter += 1;
        } else {
            self.mm_counter += 1;
        }
    }

    pub fn update_mapq_count(&mut self, rec1: &bam::Record, rec2: &bam::Record) {
        if rec1.mapq() == 0 && rec2.mapq() == 0 { // Both of alignments are not unique
            self.mq00_counter += 1;
        } else if (rec1.mapq() == 0) != (rec2.mapq() == 0) { // One of alignments unique, another one is not
            self.mq01_counter += 1;
        } else {
            self.mq11_counter += 1;
        }
    }

    pub fn update_cis_trans_count(&mut self, recs: &[PairRecord]) {
        self.pairs_counter += recs.len() as u64;
        for rec in recs {
            if rec.name1 == rec.name2 { self.intra_counter += 1; } else { self.inter_counter += 1 }
        }
    }

    pub fn dump_stats_to_file(&self, file_path: &Path) {
        let f = File::create(file_path).expect("Problem with file");
        let mut f = BufWriter::new(f);

        writeln!(f, "Total number of good aligned reads {}", self.read_counter).expect("Problem with writing file");
        writeln!(f, "Total number of good alignment records {}", self.alignment_counter).expect("Problem with writing file");
        writeln!(f, "Total number of saved HiC pairs {}", self.pairs_counter).expect("Problem with writing file");

        writeln!(f, "\nMapping Statistics").expect("Problem with writing file");
        writeln!(f, "\tBoth reads in pair are unmapped {}", self.nn_counter).expect("Problem with writing file");
        writeln!(f, "\tOne of the reads in pair is unmapped {}", self.nm_counter).expect("Problem with writing file");
        writeln!(f, "\tBoth reads in pair are mapped {}", self.mm_counter).expect("Problem with writing file");

        writeln!(f, "\nAlignment Type Statistics").expect("Problem with writing file");
        writeln!(f, "\tPrimary alignments {}", self.prime_counter).expect("Problem with writing file");
        writeln!(f, "\tSupplementary alignments {}", self.supp_counter).expect("Problem with writing file");
        writeln!(f, "\tSecondary alignments {}", self.sec_counter).expect("Problem with writing file");

        writeln!(f, "\nPair Alignment Statistics").expect("Problem with writing file");
        writeln!(f, "\tBoth reads have non-unique mapping {}", self.mq00_counter).expect("Problem with writing file");
        writeln!(f, "\tOne of the reads has non-unique mapping {}", self.mq01_counter).expect("Problem with writing file");
        writeln!(f, "\tBoth reads have unique mapping {}", self.mq11_counter).expect("Problem with writing file");

        writeln!(f, "\nHi-C Pair Statistics").expect("Problem with writing file");
        writeln!(f, "\tIntra (cis-) tig interaction {}", self.intra_counter).expect("Problem with writing file");
        writeln!(f, "\tInter (trans-) tig interaction {}", self.inter_counter).expect("Problem with writing file");
        writeln!(f, "").expect("Problem with writing file");
        writeln!(f, "\tUnique pairs {}", self.uu_pair_counter).expect("Problem with writing file");
        writeln!(f, "\tSimple pairs {}", self.uw_pair_counter).expect("Problem with writing file");
        writeln!(f, "\tvpairs {}", self.ww_pair_counter).expect("Problem with writing file");

        f.flush().expect("Problem with flushing");
    }
} 
