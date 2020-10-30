use std::error::Error;
use std::io;
use std::path::Path;

use fern;
use clap::{Arg, App, SubCommand};
use hic_convertor::{convert_bam_to_pairs, deduplicate_pairs, sort_pairs};


fn setup_logging(verbosity: u64, log_file: &Path) -> Result<(), fern::InitError> {
    let mut base_config = fern::Dispatch::new();

    base_config = match verbosity {
        0 => base_config.level(log::LevelFilter::Info),
        1 => base_config.level(log::LevelFilter::Debug),
        _ => base_config.level(log::LevelFilter::Trace),
    };

    let file_config = fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                record.target(),
                record.level(),
                message
            ))
        })
        .chain(fern::log_file(log_file)?);

    let stdout_config = fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "[{}][{}][{}] {}",
                chrono::Local::now().format("%H:%M"),
                record.target(),
                record.level(),
                message
            ))
        })
        .chain(io::stdout());

    base_config
        .chain(file_config)
        .chain(stdout_config)
        .apply()?;

    Ok(())
}

fn log_level_arg() -> Arg<'static, 'static> {
    Arg::<'static, 'static>::with_name("log_level")
        .short("l")
        .long("log_level")
        .value_name("NUM")
        .takes_value(true)
        .required(false)
        .help("Verbosity of logging (0 - 3)")
}

fn pairs_arg(hm: &'static str) -> Arg<'static, 'static> {
    Arg::<'static, 'static>::with_name("pairs")
        .short("p")
        .long("pairs")
        .value_name("FILE")
        .takes_value(true)
        .required(true)
        .help(hm)
}

fn out_pairs_arg(hm: &'static str) -> Arg<'static, 'static> {
    Arg::<'static, 'static>::with_name("out_pairs")
        .short("o")
        .long("out_pairs")
        .value_name("FILE")
        .takes_value(true)
        .required(true)
        .help(hm)
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("convertor")
        .version("0.1.0")
        .author("Pavel Avdeyev")
        .about("hic-convertor converts BAM files with Hi-C reads to Hi-C pairs. \
                It is also sorts and deduplicates obtained Hi-C reads.")
        .subcommand(
            SubCommand::with_name("convert")
                .about("Convert bam to pairs.")
                .arg(
                    Arg::with_name("bam")
                        .short("b")
                        .long("bam")
                        .value_name("FILE")
                        .takes_value(true)
                        .required(true)
                        .help("Path to alignments in bam format.")
                )
                .arg( pairs_arg("Path to file in pairs format.") )
                .arg( Arg::with_name("stats")
                    .short("s")
                    .long("stats")
                    .value_name("FILE")
                    .takes_value(true)
                    .required(true)
                    .help("Path to file with statistic.") )
                .arg(
                    Arg::with_name("graph")
                        .short("g")
                        .long("graph")
                        .value_name("FILE")
                        .takes_value(true)
                        .required(false)
                        .help("Path to graph in gfa format.")
                )
                .arg( log_level_arg() )
        )
        .subcommand(
            SubCommand::with_name("sort")
                .about("Sort pairs file using sort command (see man sort).")
                .arg( pairs_arg("Path to file with pairs.") )
                .arg( out_pairs_arg("Path to file with sorted pairs.") )
                .arg(
                    Arg::with_name("mem")
                        .short("m")
                        .long("memory")
                        .value_name("STR")
                        .takes_value(true)
                        .required(false)
                        .help("The amount of RAM memory for sorting.")
                )
                .arg(
                    Arg::with_name("nproc")
                        .short("t")
                        .long("nproc")
                        .value_name("NUM")
                        .takes_value(true)
                        .required(false)
                        .help("Number of processes for sorting.")
                )
                .arg(
                    Arg::with_name("tmpdir")
                        .short("d")
                        .long("tmpdir")
                        .value_name("PATH")
                        .takes_value(true)
                        .required(false)
                        .help("Directory for storing temporary files.")
                )
                .arg(log_level_arg() )
        )
        .subcommand(
            SubCommand::with_name("dedup")
                .about("Remove duplicated Hi-C reads from file.")
                .arg( pairs_arg("Path to file with pairs.") )
                .arg( out_pairs_arg("Path to file with deduplicated pairs.") )
                .arg(log_level_arg() )
        )
        .get_matches();

    match matches.subcommand() {
        ("convert", Some(convert_matches)) => {
            setup_logging(1, "convert.log".as_ref()).expect("failed to initialize logging.");
            let bam_file = convert_matches.value_of("bam").expect("Input bam file must be provided.");
            let pairs_file = convert_matches.value_of("pairs").expect("Output pairs file must be provided.");
            let stat_file = convert_matches.value_of("stats").expect("Output stat file must be provided.");
            match convert_matches.value_of("graph") {
                None =>  convert_bam_to_pairs(Path::new(bam_file), Path::new(pairs_file), Path::new(stat_file), None)?,
                Some(_) => convert_bam_to_pairs(Path::new(bam_file), Path::new(pairs_file), Path::new(stat_file), None)?,
            }
        },
        ("sort", Some(sort_matches)) => {
            setup_logging(1, "sort.log".as_ref()).expect("failed to initialize logging.");
            let in_file = sort_matches.value_of("pairs").expect("Input pairs file must be provided.");
            let out_file = sort_matches.value_of("out_pairs").expect("Output pairs file must be provided.");
            let nproc: u8 = sort_matches.value_of("nproc").unwrap_or("4").parse().unwrap();
            let mem= sort_matches.value_of("mem").unwrap_or("2G");
            let tmpdir = sort_matches.value_of("tmpdir");
            sort_pairs(Path::new(in_file), Path::new(out_file), nproc, mem, tmpdir)?;

        },
        ("dedup", Some(dedup_matches)) => {
            setup_logging(1, "dedup.log".as_ref()).expect("failed to initialize logging.");
            let in_file = dedup_matches.value_of("pairs").expect("Input pairs file must be provided.");
            let out_file = dedup_matches.value_of("out_pairs").expect("Output pairs file must be provided.");
            deduplicate_pairs(Path::new(in_file), Path::new(out_file));
        }
        ("", None) => eprintln!("No subcommands were provided. See help for available one."),
        _ => unreachable!(),
    };
    Ok(())
}