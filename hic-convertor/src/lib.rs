use std::path::Path;
use std::io;
use log::info;

mod pair_record;
pub mod convertor;
mod sort;
mod dedup;

pub fn convert_bam_to_pairs(bam_file: &Path, pairs_file: &Path,
                            stat_file: &Path, _graph_file: Option<&Path>) -> io::Result<()> {
    info!("Starting converting .bam to .pairs...");
    let mut converter = convertor::Converter::new(bam_file, None, pairs_file);
    converter.convert()?;
    converter.save_statistic(stat_file);
    info!("Converting .bam to .pairs is complete.");
    Ok(())
}

pub fn sort_pairs(in_file: &Path, out_file: &Path, nproc: u8, mem: &str, tmpdir: Option<&str>) -> io::Result<()> {
    info!("Starting sorting {}...", in_file.to_str().unwrap());
    sort::sort_pairs(in_file.to_str().unwrap(), out_file.to_str().unwrap(), nproc, mem, tmpdir);
    info!("Sorting results saved into {}...", out_file.to_str().unwrap());
    Ok(())
}

pub fn deduplicate_pairs(in_file: &Path, out_file: &Path) {
    info!("Starting deduplicating pairs...");
    dedup::deduplicate_pairs(in_file, out_file, 3);
    info!("Deduplicating is complete.");
}
