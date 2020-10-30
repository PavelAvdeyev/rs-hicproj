use std::process::Command;
use std::io::{self, Write};

use log::info;

use super::pair_record;

pub fn sort_pairs(pairs_path: &str, output_path: &str, nproc: u8, memory: &str, tmpdir: Option<&str>) -> io::Result<()> {
    info!("Sorting pairs in file {}", pairs_path);

    info!("Starting sorting....");

    let mut sort_c = Command::new("sort");

    sort_c.arg("-k")
            .arg(format!("{0},{0}", pair_record::COL_TIG1 + 1))
            .arg("-k")
            .arg(format!("{0},{0}", pair_record::COL_TIG2 + 1))
            .arg("-k")
            .arg(format!("{0},{0}n", pair_record::COL_POS1 + 1))
            .arg("-k")
            .arg(format!("{0},{0}n", pair_record::COL_POS2 + 1))
            .arg("--stable")
            .arg("--field-separator=\t")
            .arg(format!("--parallel={}", nproc));

    if let Some(td) = tmpdir {
        sort_c.arg(format!("--temporary-directory={}", td));
    }

    sort_c.arg("-S")
            .arg(memory)
            .arg("-o")
            .arg(output_path)
            .arg(pairs_path);

    let sort_c = sort_c.output()?;
    io::stdout().write_all(&sort_c.stdout).unwrap();
    io::stderr().write_all(&sort_c.stderr).unwrap();
    info!("Command status: {}", sort_c.status);

    info!("Done with sorting pairs.");
    
    Ok(())
}
