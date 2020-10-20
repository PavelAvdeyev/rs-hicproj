use ndarray::{s, Array1, Zip};
use std::iter::FromIterator;
use hdf5::types;
use ascii::AsciiString;
use std::path::Path;

pub type PixelT = (u32, u32, u32);

#[derive(Clone,Debug)]
pub struct MatrixReader {
    file: hdf5::File,
}

impl MatrixReader {
    pub fn new(file_path: &Path) -> hdf5::Result<MatrixReader> {
        Ok(MatrixReader {
            file: hdf5::File::open_rw(file_path)?,
        })
    }

    pub fn get_n_chroms(&self) -> hdf5::Result<usize> {
        let grp = self.file.group("chroms")?;
        Ok(grp.dataset("name")?.size())
    }

    pub fn read_resolutions(&self) -> hdf5::Result<Vec<u32>> {
        let grp = self.file.group("/resolutions/")?;
        let mmn = grp.member_names()?;
        Ok(Vec::from_iter(mmn.into_iter().map(|s_res| {
            s_res.parse::<u32>().expect("Problem with converting resolution to u32")
        })))
    }

    pub fn read_chroms_info(&self) -> hdf5::Result<(Array1<AsciiString>, Array1<u64>)> {
        let tig_orders = self.read_chrom_orders()?;
        let tig_lengths = self.read_chrom_lengths()?;
        Ok((tig_orders, tig_lengths))
    }

    pub fn get_res_group_reader(&self, res: u32) -> hdf5::Result<ResGrpReader> {
        let root = self.file.group(format!("/resolutions/{}", res).as_ref())?;
        ResGrpReader::new(root)
    }

    pub fn read_chrom_orders(&self) -> hdf5::Result<Array1<AsciiString>> {
        let grp = self.file.group("chroms")?;
        let tig_orders= read_dataset::<types::VarLenAscii>(&grp, "name")?;
        let tig_orders = Array1::from_iter(tig_orders.iter()
            .map(|x| {AsciiString::from_ascii(x.as_bytes()).unwrap()}));
        Ok(tig_orders)
    }

    pub fn read_chrom_lengths(&self) -> hdf5::Result<Array1<u64>> {
        let grp = self.file.group("chroms")?;
        read_dataset::<u64>(&grp, "length")
    }
}

#[derive(Clone,Debug)]
pub struct ResGrpReader {
    root: hdf5::Group,
}

impl ResGrpReader {
    fn new(grp: hdf5::Group) -> hdf5::Result<ResGrpReader> {
        Ok(ResGrpReader {
            root: grp,
        })
    }

    pub fn get_root(&self) -> &hdf5::Group {
        &self.root
    }

    pub fn get_n_bins(&self) -> hdf5::Result<usize> {
        let grp = self.root.group("bins")?;
        Ok(grp.dataset("chrom")?.size())
    }

    pub fn get_n_pixels(&self) -> hdf5::Result<usize> {
        let grp = self.root.group("pixels")?;
        Ok(grp.dataset("bin1_id")?.size())
    }

    pub fn read_indices(&self) -> hdf5::Result<(Array1<u32>, Array1<u32>)> {
        let tig_offsets = self.read_chrom_offsets()?;
        let bin_offsets = self.read_bin_offsets()?;
        Ok((tig_offsets, bin_offsets))
    }

    pub fn read_bin_coords(&self) -> hdf5::Result<Array1<(u32, u32)>> {
        let chroms = self.read_bin_table_chr_ids()?;
        let starts = self.read_bin_table_starts()?;
        Ok(Zip::from(&chroms).and(&starts).apply_collect(|i, s| {(*i, *s)}))
    }

    pub fn read_pixels(&self) -> hdf5::Result<(Array1<u32>, Array1<u32>, Array1<u32>)> {
        let bin1 = self.read_pixels_bin1()?;
        let bin2 = self.read_pixels_bin2()?;
        let count = self.read_pixels_count()?;
        Ok((bin1, bin2, count))
    }

    pub fn read_pixel_chunk(&self, start: usize, end: usize) -> hdf5::Result<(Array1<u32>, Array1<u32>, Array1<u32>)> {
        let bin1 = self.read_pixels_slice_bin1(start, end)?;
        let bin2 = self.read_pixels_slice_bin2(start, end)?;
        let count = self.read_pixels_slice_count(start, end)?;
        Ok((bin1, bin2, count))
    }

    // pub fn read_pixels(&self) -> hdf5::Result<Array1<PixelT>> {
    //     let bin1 = self.read_pixels_bin1()?;
    //     let bin2 = self.read_pixels_bin2()?;
    //     let count = self.read_pixels_count()?;
    //     Ok(Zip::from(&bin1).and(&bin2).and(&count).apply_collect(|b1, b2, c| { (*b1, *b2, *c) }))
    // }

    pub fn read_chrom_offsets(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("indexes")?;
        let tig_offsets = read_dataset::<u32>(&grp, "chrom_offset")?;
        Ok(tig_offsets)
    }

    pub fn read_bin_offsets(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("indexes")?;
        let bin_offsets = read_dataset::<u32>(&grp, "bin1_offset")?;
        Ok(bin_offsets)
    }

    pub fn read_bin_table_chr_ids(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("bins")?;
        read_dataset::<u32>(&grp, "chrom")
    }

    pub fn read_bin_table_starts(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("bins")?;
        read_dataset::<u32>(&grp, "start")
    }

    pub fn read_bin_table_ends(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("bins")?;
        read_dataset::<u32>(&grp, "end")
    }

    pub fn read_bin_table_weights(&self) -> hdf5::Result<Array1<f64>> {
        let grp = self.root.group("bins")?;
        read_dataset::<f64>(&grp, "weight")
    }

    pub fn read_pixels_bin1(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("pixels")?;
        read_dataset::<u32>(&grp, "bin1_id")
    }

    pub fn read_pixels_slice_bin1(&self, start: usize, end: usize) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("pixels")?;
        read_dataset_slice::<u32>(&grp, "bin1_id", start, end)
    }

    pub fn read_pixels_bin2(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("pixels")?;
        read_dataset::<u32>(&grp, "bin2_id")
    }

    pub fn read_pixels_slice_bin2(&self, start: usize, end: usize) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("pixels")?;
        read_dataset_slice::<u32>(&grp, "bin2_id", start, end)
    }

    pub fn read_pixels_count(&self) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("pixels")?;
        read_dataset::<u32>(&grp, "count")
    }

    pub fn read_pixels_slice_count(&self, start: usize, end: usize) -> hdf5::Result<Array1<u32>> {
        let grp = self.root.group("pixels")?;
        read_dataset_slice::<u32>(&grp, "count", start, end)
    }
}

pub fn read_dataset_slice<T: hdf5::H5Type>(grp: &hdf5::Group, name: &str, start: usize, end: usize)
                         -> hdf5::Result<Array1<T>> {
    let dts = grp.dataset(name)?;
    dts.read_slice_1d(s![start..end])
}

pub fn read_dataset<T: hdf5::H5Type>(grp: &hdf5::Group, name: &str) -> hdf5::Result<Array1<T>> {
    let dts = grp.dataset(name)?;
    dts.read_1d::<T>()
}

