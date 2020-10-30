use std::error::Error;
use std::io;
use std::path::Path;

use fern;
use clap::{Arg, App, SubCommand};
use hic_matrix::{zoom, Strategy, balance, create_matrix_from_pairs};


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

fn matrix_arg() -> Arg<'static, 'static> {
    Arg::<'static, 'static>::with_name("matrix")
        .short("m")
        .long("matrix")
        .value_name("FILE")
        .takes_value(true)
        .required(true)
        .help("Matrix file in specific hdf5 format.")
}

fn rslns_arg(h: &'static str) -> Arg<'static, 'static> {
    Arg::<'static, 'static>::with_name("rslns")
        .short("r")
        .long("rslns")
        .multiple(true)
        .use_delimiter(true)
        .value_terminator(";")
        .takes_value(true)
        .value_name("INT")
        .required(true)
        .help(h)
}

fn strategy_arg() -> Arg<'static, 'static> {
    Arg::<'static, 'static>::with_name("strategy")
        .short("s")
        .long("strategy")
        .possible_values(&["ICGW", "LEN"])
        .takes_value(true)
        .required(false)
        .help("Balancing strategy:. ICGW - iterative correction genome-wide, LEN - resolution size")
}

fn parse_rslns_arg(arg: Option<clap::Values>) -> Vec<u32> {
    arg.expect("List of resolutions must be provided")
        .into_iter()
        .map(|x| { x.parse().unwrap() })
        .collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("hic-matrix")
        .version("0.1.0")
        .author("Pavel Avdeyev")
        .about("The minimum code for creating/balancing/zooming Hi-C matrices.")
        .subcommand(
            SubCommand::with_name("build")
                .arg(
                    Arg::with_name("pairs")
                        .short("p")
                        .long("pairs")
                        .value_name("FILE")
                        .takes_value(true)
                        .required(true)
                        .help("A file with Hi-C pairs. File must be tab separated.
                                    1 col - read name,
                                    2 col - first contig,
                                    3 col - first coordinate,
                                    4 col - second contig,
                                    5 col - second coordinate,
                                    6 col - first strand,
                                    7 col - second strand.")
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
                .arg( rslns_arg("List of matrix resolutions") )
                .arg( matrix_arg() )
                .arg( strategy_arg() )
        ).subcommand(
            SubCommand::with_name("balance")
                .arg( matrix_arg() )
                .arg( rslns_arg("List of resolutions for balancing (must exist).") )
                .arg( strategy_arg() )
        )
        .subcommand(
            SubCommand::with_name("zoom")
                .arg(matrix_arg() )
                .arg( rslns_arg("New matrix resolutions (it must be divisable by existed resolutions in matrix)") )
        )
        .get_matches();


    match matches.subcommand() {
        ("build", Some(build_matches)) => {
            setup_logging(1, "matrix.log".as_ref()).expect("failed to initialize logging.");
            let pairs_file = Path::new(build_matches.value_of("pairs").expect("Pairs file must be provided."));
            let tig_length_file = Path::new(build_matches.value_of("lengths").expect("File with contig lengths must be provided."));
            let rslns: Vec<u32> = parse_rslns_arg(build_matches.values_of("rslns") );
            let matrix_file = Path::new(build_matches.value_of("matrix").expect("Matrix file must be provided."));
            let strategy = Strategy::from_option(build_matches.value_of("strategy"));
            create_matrix_from_pairs(pairs_file, tig_length_file, matrix_file, &rslns, &strategy)?;
        }
        ("balance", Some(bal_matches)) => {
            setup_logging(1, "matrix.log".as_ref()).expect("failed to initialize logging.");
            let matrix_file = Path::new(bal_matches.value_of("matrix").expect("Matrix file must be provided."));
            let rslns: Vec<u32> = parse_rslns_arg(bal_matches.values_of("rslns") );
            let strategy = Strategy::from_option(bal_matches.value_of("strategy"));
            balance(matrix_file, &rslns, &strategy)?;
        }
        ("zoom", Some(zoom_matches)) => {
            setup_logging(1, "matrix.log".as_ref()).expect("failed to initialize logging.");
            let matrix_file = Path::new(zoom_matches.value_of("matrix").expect("Matrix file must be provided."));
            let rslns: Vec<u32> = parse_rslns_arg(zoom_matches.values_of("rslns") );
            zoom(matrix_file, &rslns)?;
        }
        ("", None) => eprintln!("None subcommand was used. See help for available one."),
        _ => unreachable!(),
    }
    Ok(())
}
