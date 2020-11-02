use std::error::Error;
use std::iter::FromIterator;
use std::path::Path;


use ascii::AsciiString;
use ndarray::{Array1, ArrayView1};
use hdf5::types;

use super::reader::PixelT;
use super::builders::pair_builder::PairsBuilder;
use super::builders::res_grp_builder::ResGrpBuilder;

enum MatrixWriterMode {
    Write,
    Append
}

pub struct MatrixWriter {
    wrt_mode: MatrixWriterMode,
    file: hdf5::File,
}

impl MatrixWriter {
    pub fn new_in_writing_mode(filename: &Path) -> hdf5::Result<MatrixWriter> {
        MatrixWriter::new(filename, MatrixWriterMode::Write)
    }

    pub fn new_in_appending_mode(filename: &Path) -> hdf5::Result<MatrixWriter> {
        MatrixWriter::new(filename, MatrixWriterMode::Append)
    }

    pub fn get_file_handler(&self) -> &hdf5::File {
        &self.file
    }

    pub fn write_matrix(&self, builder: &PairsBuilder) -> Result<(), Box<dyn Error>> {
        self.write_chroms_group(builder.tig_names_view(), builder.tig_lengths_view())?;
        self.write_resolution_group(builder)?;
        Ok(())
    }

    pub fn write_resolution_group(&self, builder: &impl ResGrpBuilder) -> Result<(), Box<dyn Error>> {
        let grp = self.file.create_group(format!("resolutions/{}", builder.get_resolution()).as_ref())?;
        ResGrpWriter::write_resolution_group(&grp, builder)?;
        Ok(())
    }

    pub fn write_balancing_weights(&self, res: u32, weights: ArrayView1<f64>) -> hdf5::Result<()> {
        match self.wrt_mode {
            MatrixWriterMode::Write => {
                return Err(hdf5::Error::Internal(String::from("File opened in non-appending mode")));
            }
            MatrixWriterMode::Append => {
                let root = self.file.group(format!("resolutions/{}", res).as_ref())?;
                ResGrpWriter::write_balancing_weights(&root, weights)?;
            }
        };
        Ok(())
    }

    fn new(filename: &Path, wrt_mode: MatrixWriterMode) -> hdf5::Result<MatrixWriter> {
        match wrt_mode {
            MatrixWriterMode::Write => Ok(MatrixWriter {
                file: hdf5::File::create(filename)?,
                wrt_mode
            }),
            MatrixWriterMode::Append => Ok(MatrixWriter {
                file: hdf5::File::open_rw(filename)?,
                wrt_mode
            }),
        }
    }

    fn write_chroms_group(&self, tig_order: ArrayView1<AsciiString>, tig_lengths: ArrayView1<u64>) -> hdf5::Result<()> {
        let grp = self.file.create_group("chroms")?;

        let tig_orders = Array1::from_iter(
            tig_order.iter()
                .map(|x| {types::VarLenAscii::from_ascii(x.as_bytes()).unwrap()} )
        );
        write_dataset(&grp, "name", tig_orders.len(), tig_orders.view())?;
        write_dataset(&grp, "length", tig_lengths.len(), tig_lengths)?;

        Ok(())
    }
}

struct ResGrpWriter {}

impl ResGrpWriter {

    fn write_balancing_weights(grp: &hdf5::Group, weights: ArrayView1<f64>) -> hdf5::Result<()> {
        let grp = grp.group("bins")?;
        match grp.dataset("weight") {
            Ok(dts) => {
                dts.resize(weights.len())?;
                dts.write(weights);
            }
            _ => write_dataset(&grp, "weight", weights.len(), weights)?
        };

        Ok(())
    }

    fn write_resolution_group(grp: &hdf5::Group, builder: &impl ResGrpBuilder) -> Result<(), Box<dyn Error>> {
        // Writing indexes
        let pixels = ResGrpWriter::write_index_group(grp, builder)?;

        // Saving pixels
        ResGrpWriter::consume_and_write_pixels(grp, pixels)?;

        // Saving bin information
        ResGrpWriter::write_bins_description(grp, builder)?;

        Ok(())
    }

    fn write_bins_description(grp: &hdf5::Group, builder: &impl ResGrpBuilder) -> hdf5::Result<()> {
        let grp = grp.create_group("bins")?;
        let (chrs, starts, ends) = builder.get_bin_table();
        write_dataset(&grp, "chrom", chrs.len(), chrs.view())?;
        write_dataset(&grp, "start", starts.len(), starts.view())?;
        write_dataset(&grp, "end", ends.len(), ends.view())?;
        Ok(())
    }

    fn write_index_group(grp: &hdf5::Group, builder: &impl ResGrpBuilder) -> Result<Vec<PixelT>, Box<dyn Error>> {
        let grp = grp.create_group("indexes")?;

        let tig_ofssets = builder.get_tig_offsets_view();
        let pixels = builder.get_pixels()?;
        write_dataset(&grp,"chrom_offset",tig_ofssets.len(), tig_ofssets)?;
        let bin_offsets = builder.get_bin_offsets(&pixels);
        write_dataset(&grp,"bin1_offset",bin_offsets.len(), bin_offsets.view())?;

        Ok(pixels)
    }

    fn consume_and_write_pixels(grp: &hdf5::Group, pixels: Vec<PixelT>) -> hdf5::Result<()> {
        let grp = grp.create_group("pixels")?;

        let mut bin1_ids: Array1<u32> = Array1::default(pixels.len());
        let mut bin2_ids: Array1<u32> = Array1::default(pixels.len());
        let mut counts: Array1<u32> = Array1::default(pixels.len());
        pixels.into_iter().enumerate().for_each(|(i, info)| {
            let (bin1_id, bin2_id, count) = info;
            bin1_ids[i] = bin1_id;
            bin2_ids[i] = bin2_id;
            counts[i] = count;
        });

        write_dataset(&grp,"bin1_id",bin1_ids.len(), bin1_ids.view())?;
        write_dataset(&grp,"bin2_id",bin2_ids.len(), bin2_ids.view())?;
        write_dataset(&grp,"count",counts.len(), counts.view())?;

        Ok(())
    }
}

pub fn write_dataset<Q: hdf5::H5Type>(grp: &hdf5::Group, name: &str, shape: usize, ar: ArrayView1<Q>)
    -> hdf5::Result<()> {
    let dts = grp.new_dataset::<Q>().create(name, shape)?;
    dts.write(ar)?;
    Ok(())
}


// pub struct ResGrpWriter<'a, T: ResGrpBuilder> {
//     builder: &'a T,
//     root: hdf5::Group
// }

// impl<'a> ResGrpWriter<'a, PairsBuilder> {
//     pub fn from_builder(builder: &'a PairsBuilder) -> hdf5::Result<ResGrpWriter<'a, PairsBuilder>> {
//         let file = hdf5::File::create(builder.get_matrix_file())?;
//         Ok(ResGrpWriter {
//             builder,
//             root: file.create_group(format!("resolutions/{}", builder.get_resolution()).as_ref())?,
//         })
//     }
// }
//
// impl<'a> ResGrpWriter<'a, ZoomBuilder<'a>> {
//     pub fn from_zoomifier(builder: &'a ZoomBuilder<'a>) -> hdf5::Result<ResGrpWriter<'a, ZoomBuilder<'a>>> {
//         let file = hdf5::File::open_rw(builder.get_matrix_file())?;
//         Ok(ResGrpWriter {
//             builder,
//             root: file.create_group(format!("resolutions/{}", builder.get_resolution()).as_ref())?,
//         })
//     }
// }


// impl<'a, T: ResGrpBuilder> ResGrpWriter<'a, T> {
//     pub fn write(&self) -> Result<(), Box<dyn Error>> {
//         self.write_from_builder()
//     }
//
//     pub fn write_weights_for_balancing(matrix_file: &'a Path, resolution: u32, weights: &Array1<f64>)
//         -> Result<(), Box<dyn Error>> {
//         let file = hdf5::File::open_rw(matrix_file)?;
//         let root = file.group(format!("resolutions/{}", resolution).as_ref())?;
//         let grp = root.group("bins")?;
//         write_dataset(&grp, "weight", weights.len(), weights.view())?;
//         Ok(())
//     }
//
//     fn write_from_builder(&self) -> Result<(), Box<dyn Error>> {
//         // Writing indexes
//         let pixels = self.write_index_group()?;
//
//         // Saving pixels
//         self.consume_and_write_pixels(pixels)?;
//
//         // Saving bin information
//         self.write_bins_description()?;
//
//         // saving tig lengths
//         //self.write_chroms_group()?;
//
//         // Writing info
//         // self.write_meta_info()?;
//
//         Ok(())
//     }
//
//     fn write_bins_description(&self) -> hdf5::Result<()> {
//         let grp = self.root.create_group("bins")?;
//         let (chrs, starts, ends) = self.builder.get_bin_table();
//         write_dataset(&grp, "chrom", chrs.len(), chrs.view())?;
//         write_dataset(&grp, "start", starts.len(), starts.view())?;
//         write_dataset(&grp, "end", ends.len(), ends.view())?;
//         Ok(())
//     }
//
//     fn write_index_group(&self) -> Result<Vec<PixelT>, Box<dyn Error>> {
//         let grp = self.root.create_group("indexes")?;
//
//         let tig_ofssets = self.builder.get_tig_offsets_view();
//         let pixels = self.builder.build_pixels()?;
//         write_dataset(&grp,"chrom_offset",tig_ofssets.len(), tig_ofssets)?;
//         let bin_offsets = self.builder.build_bin_offsets(&pixels);
//         write_dataset(&grp,"bin1_offset",bin_offsets.len(), bin_offsets.view())?;
//
//         Ok(pixels)
//     }
//
//     fn consume_and_write_pixels(&self, pixels: Vec<PixelT>) -> hdf5::Result<()> {
//         let grp = self.root.create_group("pixels")?;
//
//         let mut bin1_ids: Array1<u32> = Array1::default(pixels.len());
//         let mut bin2_ids: Array1<u32> = Array1::default(pixels.len());
//         let mut counts: Array1<u32> = Array1::default(pixels.len());
//         pixels.into_iter().enumerate().for_each(|(i, info)| {
//             let (bin1_id, bin2_id, count) = info;
//             bin1_ids[i] = bin1_id;
//             bin2_ids[i] = bin2_id;
//             counts[i] = count;
//         });
//
//         write_dataset(&grp,"bin1_id",bin1_ids.len(), bin1_ids.view())?;
//         write_dataset(&grp,"bin2_id",bin2_ids.len(), bin2_ids.view())?;
//         write_dataset(&grp,"count",counts.len(), counts.view())?;
//
//         Ok(())
//     }
// }
