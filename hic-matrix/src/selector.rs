use ndarray::{s, Array1};
use super::{utils, reader::ResGrpReader};
use std::{mem, iter};
use itertools::{Itertools, izip};
use std::iter::FromIterator;

#[derive(Clone,Debug)]
pub struct Selector2D {
    bin_offsets: Array1<u32>,
    biases: Array1<f64>,
    reader: ResGrpReader
}

impl Selector2D {
    pub fn new(reader: ResGrpReader) -> hdf5::Result<Selector2D> {
        Ok(Selector2D {
            bin_offsets: reader.read_bin_offsets()?,
            biases: reader.read_bin_table_weights()?,
            reader
        })
    }

    pub fn get_balanced_submatrix(&self, i0: usize, i1: usize, j0: usize, j1: usize)
        -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<f64>)> {
        let (is, js, vs) = self.get_rectangle(i0, i1, j0, j1)?;
        let bvs= Vec::from_iter(izip!(is.iter(), js.iter(), vs.iter()).map(|(&b1, &b2, &v)| {
            (v  as f64) * self.biases[b1 as usize] * self.biases[b2 as usize]
        }));
        Ok((is, js, bvs))
    }

    pub fn get_raw_submatrix(&self, i0: usize, i1: usize, j0: usize, j1: usize)
        -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<u32>)> {
        self.get_rectangle(i0, i1, j0, j1)
    }

    fn get_rectangle(&self, mut i0: usize, mut i1: usize, mut j0: usize, mut j1: usize)
        -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<u32>)> {
        let mut is;
        let mut js;
        let vs;

        if (i0, i1) == (j0, j1) {
            // symmetric query
            let (tis, tjs, tvs) = self.get_triu_nnz_bins(i0, i1, i0, i1)?;
            let nodiags = Vec::from_iter(tis.iter().zip(tjs.iter()).map(|(&a, &b)| a != b));

            let tis_nd = utils::get_vec_wrt_predicate(&nodiags, &tis);
            let tjs_nd = utils::get_vec_wrt_predicate(&nodiags, &tjs);
            let tvs_nd = utils::get_vec_wrt_predicate(&nodiags, &tvs);

            is = [&tis[..], &tjs_nd[..]].concat(); // np.r_[i, j[nodiag]]
            js = [&tjs[..], &tis_nd[..]].concat(); // np.r_[j, i[nodiag]]
            vs = [&tvs[..], &tvs_nd[..]].concat(); // np.r_[v, v[nodiag]]
        } else {
            // asymmetric query
            let mut transpose = false;

            if j0 < i0 || (i0 == j0 && i1 < j1) {
                // println!("Swap occur {} {} {} {}", i0, i1, j0, j1);
                mem::swap(&mut i0, &mut j0);
                mem::swap(&mut i1, &mut j1);
                transpose = true;
            }

            // println!("{} {} {} {}", i0, i1, j0, j1);
            if !Selector2D::is_overlap(i0, i1, j0, j1) {
                // non-overlapping
                // println!("Non overlapping {} {} {} {}", i0, i1, j0, j1);
                let (tis, tjs, tvs) = self.get_triu_nnz_bins(i0, i1, j0, j1)?;
                is = tis;
                js = tjs;
                vs = tvs;
            } else if Selector2D::is_nested(i0, i1, j0, j1) {
                // nested
                // println!("Nested {} {} {} {}", i0, i1, j0, j1);
                let (ix, jx, vx) = self.get_triu_nnz_bins(i0, j0, j0, j1)?;
                let (mut jy, mut iy, mut vy) = self.get_triu_nnz_bins(j0, j1, j0, j1)?;
                let (jz, iz, vz) = self.get_triu_nnz_bins(j0, j1,j1, i1)?;

                // glue without duplication
                let nodiags = Vec::from_iter(iy.iter().zip(jy.iter()).map(|(&a, &b)| a != b));
                let iy_nd = utils::get_vec_wrt_predicate(&nodiags, &iy);
                let jy_nd = utils::get_vec_wrt_predicate(&nodiags, &jy);
                let vy_nd = utils::get_vec_wrt_predicate(&nodiags, &vy);

                iy = [&iy[..], &jy_nd[..]].concat(); //np.r_[iy, jy[nodiag]],
                jy = [&jy[..], &iy_nd[..]].concat(); //np.r_[jy, iy[nodiag]],
                vy = [&vy[..], &vy_nd[..]].concat(); //np.r_[vy, vy[nodiag]],

                is = [&ix[..], &iy[..], &iz[..]].concat(); // np.r_[ix, iy, iz]
                js = [&jx[..], &jy[..], &jz[..]].concat(); // np.r_[jx, jy, jz]
                vs = [&vx[..], &vy[..], &vz[..]].concat(); // np.r_[vx, vy, vz]
            } else if Selector2D::is_coming_before(i0, i1, j0, j1) {
                // intervals should be sorted in acceding order
                let (ix, jx, vx) = self.get_triu_nnz_bins(i0, j0, j0, i1)?;
                let (mut iy, mut jy, mut vy) = self.get_triu_nnz_bins(j0, i1, j0, i1)?;
                let (iz, jz, vz) = self.get_triu_nnz_bins(i0, i1,i1, j1)?;

                let nodiags = Vec::from_iter(iy.iter().zip(jy.iter()).map(|(&a, &b)| a != b));
                let iy_nd = utils::get_vec_wrt_predicate(&nodiags, &iy);
                let jy_nd = utils::get_vec_wrt_predicate(&nodiags, &jy);
                let vy_nd = utils::get_vec_wrt_predicate(&nodiags, &vy);

                iy = [&iy[..], &jy_nd[..]].concat(); //np.r_[iy, jy[nodiag]],
                jy = [&jy[..], &iy_nd[..]].concat(); //np.r_[jy, iy[nodiag]],
                vy = [&vy[..], &vy_nd[..]].concat(); //np.r_[vy, vy[nodiag]],

                is = [&ix[..], &iy[..], &iz[..]].concat(); // np.r_[ix, iy, iz]
                js = [&jx[..], &jy[..], &jz[..]].concat(); // np.r_[jx, jy, jz]
                vs = [&vx[..], &vy[..], &vz[..]].concat(); // np.r_[vx, vy, vz]
            } else {
                panic!("It should not happen. Please report to git.")
            }

            if transpose {
                mem::swap(&mut is, &mut js);
            }
        }

        assert!((is.len() == js.len()) && (js.len() == vs.len()));

        Ok((is, js, vs))
    }

    fn get_triu_nnz_bins(&self, i0: usize, i1: usize, j0: usize, j1: usize)
        -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<u32>)> {

        if i0 >= i1 || j0 >= j1 {
            return Ok((Vec::new(), Vec::new(), Vec::new()));
        }

        let mut is = Vec::new();
        let mut js = Vec::new();
        let mut vs = Vec::new();

        let intervals = self.bin_offsets.slice(s!(i0..=i1));
        let p0 = intervals[0] as usize;
        let p1 = intervals[intervals.len() - 1] as usize;
        let (_, bin2ids, counts) = self.reader.read_pixel_chunk(p0, p1)?;

        for (row_id, (&lo, &hi)) in (i0..i1).zip(intervals.iter().tuple_windows()) {
            let (lo, hi) = ((lo as usize) - p0, (hi as usize) - p0);

            let cur_bins = bin2ids.slice(s![lo..hi]);
            let cur_counts = counts.slice(s![lo..hi]);
            let mask = cur_bins.mapv(|x| (x >= j0 as u32) && (x < j1 as u32));
            let cols = utils::get_array_wrt_predicate(mask.view(), cur_bins);

            is.extend(iter::repeat(row_id as u32).take(cols.len()));
            js.extend(cols);
            vs.extend(utils::get_array_wrt_predicate(mask.view(), cur_counts));
        }

        assert!((is.len() == js.len()) && (js.len() == vs.len()));
        Ok((is, js, vs))
    }

    fn is_overlap(i0: usize, i1: usize, j0: usize, j1: usize) -> bool {
        (i0 <= j1) && (j0 <= i1)
    }

    fn is_nested(i0: usize, i1: usize, j0: usize, j1: usize) -> bool {
        (i0 <= j0 && j1 <= i1) || (j0 <= i0 && i1 <= j1)
    }

    fn is_coming_before(i0: usize, i1: usize, j0: usize, j1: usize) -> bool {
        (i0 < j0) && (i1 <= j1)
    }
}

// fn debug(bias: ArrayView1<i64>) {
//     use std::fs::File;
//     use std::io::{BufWriter, Write};
//     let f = File::create("test.txt").expect("Problem with file");
//     let mut f = BufWriter::new(f);
//     for x in bias.view() {
//         writeln!(f, "{}", *x).expect("Problem with writing file");
//     }
//     f.flush().expect("Problem with flushing");
// }