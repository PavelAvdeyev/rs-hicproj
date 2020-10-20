use std::path::Path;
use std::io;
use log::info;

mod pair_record;
pub mod convertor;
mod sort;
mod dedup;

pub fn convert_bam_to_pairs(bam_file: &Path, _graph_file: Option<&Path>,
                            pairs_file: &Path, stat_file: &Path) -> io::Result<()> {
    info!("Starting converting .bam to .pairs...");
    let mut converter = convertor::Converter::new(bam_file, None, pairs_file);
    converter.convert()?;
    converter.save_statistic(stat_file);
    info!("Converting .bam to .pairs is complete.");
    Ok(())
}

pub fn sort_pairs(in_file: &Path, out_file: &Path, tmp_dir: Option<&Path>, nproc: u8) -> io::Result<()> {
    info!("Starting sorting {}...", in_file.to_str().unwrap());
    match tmp_dir {
        Some(p) => sort::sort_pairs(in_file.to_str().unwrap(),
                                    out_file.to_str().unwrap(),
                                    p.to_str().unwrap(),
                                    nproc, "2G")?,
        None => sort::sort_pairs(in_file.to_str().unwrap(),
                                 out_file.to_str().unwrap(),
                                 " ",
                                 nproc, "2G")?,
    }
    info!("Sorting results saved into {}...", out_file.to_str().unwrap());
    Ok(())
}

pub fn deduplicate_pairs(in_file: &Path, out_file: &Path) {
    info!("Starting deduplicating pairs...");
    dedup::deduplicate_pairs(in_file, out_file, 3);
    info!("Deduplicating is complete.");
}

pub fn full_pipeline(bam_file: &Path, graph_file: Option<&Path>, output_dir: &Path, nproc: u8) -> io::Result<()> {
    if !output_dir.is_dir() { return Err(io::Error::new(io::ErrorKind::InvalidData, "The input is not directory")); }

    info!("Preparing pairs (convert/sort/deduplicate)...");
    let pairs_file = output_dir.join("raw.pairs");
    let stat_file = output_dir.join("stat.txt");
    convert_bam_to_pairs(bam_file, graph_file, pairs_file.as_path(), stat_file.as_path())?;

    let sort_pairs_file = output_dir.join("sorted.pairs");
    sort_pairs(pairs_file.as_path(), sort_pairs_file.as_path(), Some(output_dir), nproc)?;

    let final_file = output_dir.join("final.pairs");
    deduplicate_pairs(sort_pairs_file.as_path(), final_file.as_path());
    info!("Done with preparing pairs.");
    Ok(())
}