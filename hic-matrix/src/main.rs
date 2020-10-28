use std::error::Error;
use std::io;
use std::path::Path;

use fern;
use clap::{Arg, App, SubCommand};
use hic_matrix::{zoom_with_balancing, Strategy, balance, create_matrix_from_pairs, create_multi_matrix_from_pairs};


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

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("hic-matrix")
        .version("0.1.0")
        .author("Pavel Avdeyev")
        .about("The minimum code for creating/balancing/zooming Hi-C matrices.")
        .subcommand(
            SubCommand::with_name("all")
                .arg(
                    Arg::with_name("pairs")
                        .short("p")
                        .long("pairs")
                        .value_name("FILE")
                        .takes_value(true)
                        .required(true)
                        .help("A file with Hi-C pairs. File must be tab seperated.\
                                1 col - read name, 2 col - first contig, 3 col - first coordinate\
                                4 col - second contig, 5 col - second coordinate, \
                                6 col - first strand, 7 col - second strand.")
                )
                .arg(
                    Arg::with_name("lengths")
                        .short("l")
                        .long("lengts")
                        .value_name("FILE")
                        .takes_value(true)
                        .required(true)
                        .help("File with contig lengths. First column is tig name, \
                                second one is length.")
                )
                .arg(
                    Arg::with_name("rslns")
                        .short("r")
                        .long("rslns")
                        .multiple(true)
                        .use_delimiter(true)
                        .value_terminator(";")
                        .takes_value(true)
                        .value_name("INT")
                        .required(true)
                        .help("Resolutions that will be created from pairs.")
                )
                .arg(
                    Arg::with_name("strategy")
                        .short("s")
                        .long("strategy")
                        .possible_values(&["ICGW", "LEN"])
                        .takes_value(true)
                        .required(false)
                        .help("Balancing strategy:. ICGW - iterative correction genome-wide, LEN - resolution size")
                )
                .arg(
                    Arg::with_name("matrix")
                        .short("m")
                        .long("matrix")
                        .value_name("FILE")
                        .takes_value(true)
                        .required(true)
                        .help("Output file, where matrix will be stored in hdf5 format.")
                )
        )
        .subcommand(
            SubCommand::with_name("convert")
                        .arg(
                            Arg::with_name("pairs")
                                .short("p")
                                .long("pairs")
                                .value_name("FILE")
                                .takes_value(true)
                                .required(true)
                                .help("A file with Hi-C pairs. File must be tab seperated.\
                                1 col - read name, 2 col - first contig, 3 col - first coordinate\
                                4 col - second contig, 5 col - second coordinate, \
                                6 col - first strand, 7 col - second strand.")
                        )
                        .arg(
                            Arg::with_name("lengths")
                                .short("l")
                                .long("lengts")
                                .value_name("FILE")
                                .takes_value(true)
                                .required(true)
                                .help("File with contig lengths. First column is tig name, \
                                second one is length.")
                        )
                        .arg(
                            Arg::with_name("rsln")
                                .short("r")
                                .long("rsln")
                                .value_name("INT")
                                .takes_value(true)
                                .required(true)
                                .help("A resolution that will be created from pairs.")
                        )
                        .arg(
                            Arg::with_name("matrix")
                                .short("m")
                                .long("matrix")
                                .value_name("FILE")
                                .takes_value(true)
                                .required(true)
                                .help("Output file, where matrix will be stored in hdf5 format.")
                        )
                        .arg(
                            Arg::with_name("strategy")
                                .short("s")
                                .long("strategy")
                                .possible_values(&["ICGW", "LEN"])
                                .takes_value(true)
                                .required(false)
                                .help("Balancing strategy:. ICGW - iterative correction genome-wide, LEN - resolution size")
                        )
        )
        .subcommand(
            SubCommand::with_name("balance")
                        .arg(
                            Arg::with_name("matrix")
                                .short("m")
                                .long("matrix")
                                .value_name("FILE")
                                .takes_value(true)
                                .required(true)
                                .help("Matrix file in specific hdf5 format.")
                        )
                        .arg(
                            Arg::with_name("rsln")
                                .short("r")
                                .long("rsln")
                                .value_name("INT")
                                .takes_value(true)
                                .required(true)
                                .help("A resolution that will be balanced. It must exist.")
                        )
                        .arg(
                            Arg::with_name("strategy")
                                .short("s")
                                .long("strategy")
                                .possible_values(&["ICGW", "LEN"])
                                .takes_value(true)
                                .required(false)
                                .help("Balancing strategy:. ICGW - iterative correction genome-wide, LEN - resolution size")
                        )
        )
        .subcommand(
            SubCommand::with_name("zoom")
                .arg(
                    Arg::with_name("matrix")
                            .short("m")
                            .long("matrix")
                            .value_name("FILE")
                            .takes_value(true)
                            .required(true)
                            .help("Matrix file in specific hdf5 format.")
                )
                .arg(
                    Arg::with_name("rsln")
                        .short("r")
                        .long("rsln")
                        .value_name("INT")
                        .takes_value(true)
                        .required(true)
                        .help("New resolution, which will be created from existence ones. \
                                    It should be divisable by resolution already existing in the matrix")
                )
                .arg(
                    Arg::with_name("strategy")
                        .short("s")
                        .long("strategy")
                        .possible_values(&["ICGW", "LEN"])
                        .takes_value(true)
                        .required(false)
                        .help("Balancing strategy:. ICGW - iterative correction genome-wide, LEN - resolution size")
                )
        )
        .get_matches();

    match matches.subcommand() {
        ("all", Some(all_matches)) => {
            setup_logging(1, "matrix.log".as_ref()).expect("failed to initialize logging.");
            let pairs_file = Path::new(all_matches.value_of("pairs").unwrap());
            let tig_length_file = Path::new(all_matches.value_of("lengths").unwrap());
            let resolutions: Vec<u32> = all_matches.values_of("rslns").unwrap().into_iter().map(|x| { x.parse().unwrap() }).collect();
            let matrix_file = Path::new(all_matches.value_of("matrix").unwrap());
            let strategy = match all_matches.value_of("strategy") {
                Some(strategy) => { Strategy::from_string(strategy) },
                None => Strategy::None,
            };

            create_multi_matrix_from_pairs(pairs_file, tig_length_file, matrix_file, &resolutions, &strategy)?;
        },
        ("convert", Some(convert_matches)) => {
            setup_logging(1, "matrix.log".as_ref()).expect("failed to initialize logging.");
            let pairs_file = Path::new(convert_matches.value_of("pairs").unwrap());
            let tig_length_file = Path::new(convert_matches.value_of("lengths").unwrap());
            let rsln: u32 = convert_matches.value_of("rsln").unwrap().parse().unwrap();
            let matrix_file = Path::new(convert_matches.value_of("matrix").unwrap());
            let strategy = match convert_matches.value_of("strategy") {
                Some(strategy) => { Strategy::from_string(strategy) },
                None => Strategy::None,
            };
            create_matrix_from_pairs(pairs_file, tig_length_file, matrix_file, rsln, &strategy)?;
        },
        ("balance", Some(bal_matches)) => {
            setup_logging(1, "matrix.log".as_ref()).expect("failed to initialize logging.");
            let rsln: u32 = bal_matches.value_of("rsln").unwrap().parse().unwrap();
            let matrix_file = Path::new(bal_matches.value_of("matrix").unwrap());
            let strategy = match bal_matches.value_of("strategy") {
                Some(strategy) => { Strategy::from_string(strategy) },
                None => Strategy::None,
            };
            balance(matrix_file, rsln, &strategy)?;
        },
        ("zoom", Some(zoom_matches)) => {
            setup_logging(1, "matrix.log".as_ref()).expect("failed to initialize logging.");
            let rsln: u32 = zoom_matches.value_of("rsln").unwrap().parse().unwrap();
            let matrix_file = Path::new(zoom_matches.value_of("matrix").unwrap());
            let strategy = match zoom_matches.value_of("strategy") {
                Some(strategy) => { Strategy::from_string(strategy) },
                None => Strategy::None,
            };
            zoom_with_balancing(matrix_file, &vec![rsln], &strategy)?;
        },
        ("", None) => println!("None subcommand was used. See help for available one."),
        _ => unreachable!(),
    }
    Ok(())
}
