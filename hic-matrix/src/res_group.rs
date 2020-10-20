use ndarray::{Array1, Array2, ArrayView1};
use std::error;
use itertools::izip;
use std::iter::FromIterator;

use super::selector::Selector2D;
use super::reader::ResGrpReader;
use super::errors::{MatrixIndexError, SelectorUninitError};


#[derive(Clone,Debug)]
pub struct ResGroup {
    resolution: u32,
    n_bins: usize,
    n_pixels: usize,
    reader: ResGrpReader,
    selector: Option<Selector2D>
}

// pub struct MetaMatrixInfo {
//     bin_size: u32,
//     nbins: u32,
//     sum: u32,
//     nnz: u32,
// }


impl ResGroup {
    pub fn new(resolution: u32, reader: ResGrpReader) -> hdf5::Result<ResGroup> {
        Ok(ResGroup {
            resolution,
            n_bins: reader.get_n_bins()?,
            n_pixels: reader.get_n_pixels()?,
            reader,
            selector: None
        })
    }

    pub fn init_selector(&mut self) -> hdf5::Result<()> {
        if self.selector.is_none() {
            self.selector = Some(Selector2D::new(self.reader.clone())?)
        }
        Ok(())
    }

    pub fn get_resolution(&self) -> u32 {
        self.resolution
    }

    pub fn get_n_bins(&self) -> usize {
        self.n_bins
    }

    pub fn get_n_pixels(&self) -> usize {
        self.n_pixels
    }

    pub fn get_balanced_submatrix_as_array(&self, i0: usize, i1: usize, j0: usize, j1: usize)
        -> Result<Array2<f64>, Box<dyn error::Error>> {
        if (i0 >= i1) || (j0 >= j1) || (i1 > self.n_bins) || (j1 > self.n_bins)
            || (i0 >= self.n_bins) || (j0 >= self.n_bins) {
            return Err(MatrixIndexError.into());
        }

        if self.selector.is_none() {
            return Err(SelectorUninitError.into());
        }

        let (is, js, vs) = self.get_balanced_submatrix(i0, i1, j0, j1)?;
        let mut matrix = Array2::<f64>::zeros((i1 - i0, j1 - j0));
        for (&b1, &b2, &v) in izip!(is.iter(), js.iter(), vs.iter()) {
            let (b1, b2) = (b1 as usize - i0, b2 as usize - j0);
            if v.is_finite() { matrix[[b1, b2]] = v }
        }

        Ok(matrix)
    }

    pub fn get_balanced_column_as_array(&self, col_id: usize) -> Result<Array1<f64>, Box<dyn error::Error>> {
        self.get_balanced_row_as_array(col_id)
    }

    pub fn get_balanced_row_as_array(&self, row_id: usize) -> Result<Array1<f64>, Box<dyn error::Error>> {
        let nnz_elems = self.get_balanced_row_as_nnz_elems(row_id)?;
        let mut row = Array1::<f64>::zeros(self.n_bins);
        for (b2, v) in nnz_elems.iter() {
            row[*b2 as usize] = *v;
        }
        Ok(row)
    }

    pub fn get_balanced_row_as_nnz_elems(&self, row_id: usize) -> Result<Vec<(u32, f64)>, Box<dyn error::Error>> {
        if row_id >= self.n_bins {
            return Err(MatrixIndexError.into());
        }

        if self.selector.is_none() {
            return Err(SelectorUninitError.into());
        }

        let (_, js, vs) = self.get_balanced_submatrix(row_id, row_id + 1, 0, self.n_bins)?;
        let mut res = Vec::from_iter(js.into_iter().zip(vs.into_iter()).filter(|(_, v)| v.is_finite()));
        res.sort_by_key(|x| x.0);
        Ok(res)
    }

    pub fn get_bin_coords(&self) -> hdf5::Result<Array1<(u32, u32)>> {
        self.reader.read_bin_coords()
    }

    pub fn get_bin_chr_ids(&self) -> hdf5::Result<Array1<u32>> {
        self.reader.read_bin_table_chr_ids()
    }

    pub fn get_tigs_offsets(&self) -> hdf5::Result<Array1<u32>> {
        self.reader.read_chrom_offsets()
    }

    pub fn get_raw_pixels(&self) -> hdf5::Result<(Array1<u32>, Array1<u32>, Array1<u32>)> {
        self.reader.read_pixels()
    }

    pub fn get_balanced_pixels_range(&self, start: usize, end: usize) -> hdf5::Result<(Array1<u32>, Array1<u32>, Array1<f64>)> {
        let biases = self.reader.read_bin_table_weights()?;
        let (bins1, bins2, counts) = self.reader.read_pixel_chunk(start, end)?;
        let weights = balance_counts(biases.view(), bins1.view(), bins2.view(), counts.view());
        Ok((bins1, bins2, weights))
    }

    pub fn get_raw_pixel_iter_range(&self, start: usize, end: usize, step_l: usize) -> RawPixelIterator {
        RawPixelIterator::new(&self.reader, start, end, step_l)
    }

    pub fn get_raw_pixel_iter(&self, step_l: usize) -> RawPixelIterator {
        self.get_raw_pixel_iter_range(0, self.n_pixels, step_l)
    }

    pub fn get_balanced_pixel_iter_range(&self, start: usize, end: usize, step_l: usize) -> hdf5::Result<BalancedPixelIterator> {
        BalancedPixelIterator::new(&self.reader, start, end, step_l)
    }

    pub fn get_balanced_pixel_iter(&self, step_l: usize) -> hdf5::Result<BalancedPixelIterator> {
        self.get_balanced_pixel_iter_range(0, self.n_pixels, step_l)
    }

    fn get_balanced_submatrix(&self, i0: usize, i1: usize, j0: usize, j1: usize) -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<f64>)> {
        assert!(self.selector.is_some());
        let sel = self.selector.as_ref().unwrap();
        let (is, js, vs) = sel.get_balanced_submatrix(i0, i1, j0, j1)?;
        Ok((is, js, vs))
    }

    fn get_raw_submatrix(&self, i0: usize, i1: usize, j0: usize, j1: usize) -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<u32>)> {
        assert!(self.selector.is_some());
        let sel = self.selector.as_ref().unwrap();
        let (is, js, vs) = sel.get_raw_submatrix(i0, i1, j0, j1)?;
        // debug(&is, &js, &vs);
        Ok((is, js, vs))
    }

}

pub struct RawPixelIterator<'a> {
    reader: &'a ResGrpReader,
    current: usize,
    end: usize,
    chunksize: usize
}

impl<'a> Iterator for RawPixelIterator<'a> {
    type Item = (Array1<u32>, Array1<u32>, Array1<u32>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.end {
            let start = self.current;
            self.current = self.end.min(self.current + self.chunksize);
            self.reader.read_pixel_chunk(start, self.current).ok()
        } else {
            None
        }
    }
}

impl<'a> RawPixelIterator<'a> {
    pub fn new(reader: &'a ResGrpReader, start: usize, end: usize, chunksize: usize) -> RawPixelIterator<'a> {
        assert!(start < end);
        RawPixelIterator {
            reader,
            current: start,
            end,
            chunksize
        }
    }
}

pub struct BalancedPixelIterator<'a> {
    raw_iter: RawPixelIterator<'a>,
    biases: Array1<f64>,
}

impl<'a> Iterator for BalancedPixelIterator<'a> {
    type Item = (Array1<u32>, Array1<u32>, Array1<f64>);

    fn next(&mut self) -> Option<Self::Item> {
        self.raw_iter.next().and_then(|(bins1, bins2, counts)| {
            let weights = balance_counts(self.biases.view(), bins1.view(), bins2.view(), counts.view());
            Some((bins1, bins2, weights))
        })
    }
}

impl<'a> BalancedPixelIterator<'a> {
    pub fn new(reader: &'a ResGrpReader, start: usize, end: usize, chunksize: usize) -> hdf5::Result<BalancedPixelIterator<'a>> {
        Ok(BalancedPixelIterator {
            raw_iter: RawPixelIterator::<'a>::new(reader, start, end, chunksize),
            biases: reader.read_bin_table_weights()?
        })
    }
}

pub fn balance_counts(biases: ArrayView1<f64>, bins1: ArrayView1<u32>, bins2: ArrayView1<u32>,
                      counts: ArrayView1<u32>) -> Array1<f64> {
    Array1::from_iter(counts.iter().enumerate().map(|(i, &count)| {
        (count as f64) * biases[bins1[i] as usize] * biases[bins2[i] as usize]
    }))
}


// fn debug(is: &Vec<i64>, js: &Vec<i64>, vs: &Vec<i32>) {
//     use std::fs::File;
//     use std::io::{BufWriter, Write};
//     let f = File::create("test.txt").expect("Problem with file");
//     let mut f = BufWriter::new(f);
//     for (&x, (&y, &z)) in is.iter().zip(js.iter().zip(vs.iter())) {
//         writeln!(f, "{} {} {}", x, y, z).expect("Problem with writing file");
//     }
//     f.flush().expect("Problem with flushing");
// }

// impl MetaMatrixInfo {
//     pub fn new() -> MetaMatrixInfo {
//         MetaMatrixInfo {
//             bin_type: String::from("fixed"),
//             bin_size: 0,
//             storage_mode: String::from("symmetric-upper"),
//             nchroms: 0,
//             nbins: 0,
//             sum: 0,
//             nnz: 0,
//             genome_assembly: String::from("unknown"),
//             creation_date: Local::now().to_rfc3339(),
//             generated_by: String::from("scaff"),
//             format: String::from("HDF5::Cooler"),
//             format_version: 3.to_string(),
//             format_url: String::from("https://github.com/mirnylab/cooler"),
//         }
//     }
//
//     pub fn from(bin_size: u32, n_chroms: u32, n_bins: u32, ncont: u32, nnz: u32) -> MetaMatrixInfo {
//         let mut info = MetaMatrixInfo::new();
//         info.bin_size = bin_size;
//         info.nchroms = n_chroms;
//         info.nbins = n_bins;
//         info.sum = ncont;
//         info.nnz = nnz;
//         info
//     }
//
//     pub fn write_to_hdf5_as_attrs(&self, hdf_file: hdf5::File) {
//         //TBA
//     }
// }