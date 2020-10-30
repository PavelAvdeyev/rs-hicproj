use std::ops;
use std::iter::FromIterator;
use std::path::Path;

use ndarray::{azip, s, Array1, Zip, ArrayView1};
use ndarray_stats::SummaryStatisticsExt;
use itertools::Itertools;
use num_traits::identities;

use super::res_group::ResGroup;
use super::utils;

const CHUNKSIZE: usize = 30_000_000;

pub enum Strategy {
    ICGenomeWide,
    BinLength,
    None,
}

impl Strategy {
    pub fn from_string(s: &str) -> Strategy {
        match s {
            "ICGW" => Strategy::ICGenomeWide,
            "LEN" => Strategy::BinLength,
            _ => Strategy::None
        }
    }

    pub fn from_option(s: Option<&str>) -> Strategy {
        match s {
            Some(strategy) => { Strategy::from_string(strategy) }
            None => Strategy::None,
        }
    }
}

pub struct Balancer {
    ignore_diags: u32,
    min_nnz: u32,
    n_iters: usize,
    mad_max: f64,
    var_bound: f64,
}

impl Balancer {
    pub fn new() -> Balancer {
        Balancer {
            ignore_diags: 3,
            min_nnz: 5,
            n_iters: 400,
            mad_max: 5.0,
            var_bound: 1e-5,
        }
    }

    pub fn balance_by_resolution(&self, res_group: &ResGroup) -> Array1<f64> {
        let wght = 1.0 / (res_group.get_resolution() as f64);
        let bias = Array1::from_elem((res_group.get_n_bins(),), wght);
        bias
    }

    pub fn balance_by_ic_genomewide(&self, res_group: &ResGroup) -> Option<Array1<f64>> {
        let bias = Array1::<f64>::ones((res_group.get_n_bins(),));
        let bias = self.filter_few_nnzs(res_group, bias);
        let bias = self.filter_bins_by_mad(res_group, bias);
        let bias = self.do_iterative_corrections(res_group, bias);

        if let Some(b) = &bias {
            debug(b, format!("balance1_{}.txt", res_group.get_resolution().to_string()).as_ref());
        }

        bias
    }

    fn filter_few_nnzs(&self, res_group: &ResGroup, bias: Array1<f64>) -> Array1<f64> {
        let mut res = Array1::<u32>::zeros((res_group.get_n_bins(),));
        for (bins1, bins2, counts) in res_group.get_raw_pixel_iter(CHUNKSIZE) {
            let data = self.pipe_binarize(res_group.get_n_bins(), bins1, bins2, counts);
            res += &data;
        }

        let res = res.mapv(|m| m < self.min_nnz);
        self.filter_by_predicate(res.view(), bias)
    }

    fn filter_bins_by_mad(&self, res_group: &ResGroup, mut bias: Array1<f64>) -> Array1<f64> {
        let mut res = Array1::<f64>::zeros((res_group.get_n_bins(),));
        for (bins1, bins2, counts) in res_group.get_raw_pixel_iter(CHUNKSIZE) {
            let data = self.pipe_zeroing(res_group.get_n_bins(), bins1, bins2, counts).mapv(|x| x as f64);
            res += &data;
        }

        for (lo, hi) in res_group.get_tigs_offsets().unwrap().iter().tuple_windows() {
            let (lo, hi) = (*lo as usize, *hi as usize);

            let c_marg = res.slice(s![lo..hi]);
            let nnz_elems = utils::get_array_wrt_predicate(c_marg.mapv(|x| x > 0.0).view(), c_marg.view());

            if let Some(median) = utils::median(&nnz_elems) {
                res.slice_mut(s![lo..hi]).map_inplace(|x| *x /= median);
            } else {
                res.slice_mut(s![lo..hi]).map_inplace(|x| *x = 0.0);
            }
        }

        let mut nnz_elems = utils::get_array_wrt_predicate(res.mapv(|x| x != 0.0).view(), res.view());
        nnz_elems.iter_mut().for_each(|x| { *x = x.ln() });
        let log_nnz_med = utils::median(&nnz_elems);
        let log_nnz_dev = utils::mad(nnz_elems);
        let cutoff = log_nnz_med.zip(log_nnz_dev).map(|(med, dev)| {
            (med - self.mad_max * dev).exp()
        });

        if let Some(bound) = cutoff {
            bias = self.filter_by_predicate(res.mapv(|m| m < bound).view(), bias);
        } else {
            println!("Mad correction was not performed since problems with calculations");
        }
        bias
    }

    fn do_iterative_corrections(&self, res_group: &ResGroup, mut bias: Array1<f64>) -> Option<Array1<f64>> {
        for iteration in 0..self.n_iters {
            match self.calc_mean_and_var_of_matrix(res_group, bias.view()) {
                Some(((mean, var), mut data)) => {
                    // println!("Mean {}", mean);
                    data.map_inplace(|x| if *x == 0.0 {*x = 1.0;} else {*x /= mean;} );
                    bias = Zip::from(&bias).and(&data).apply_collect(|&b, &d| {b / d}); //TODO think about nans and infinities
                    println!("variance is {} on iteration {}", var, iteration);
                    if var < self.var_bound { break; }
                },
                _ => {
                    println!("Problem with computing mean. Abort balancing.");
                    return None;
                }
            };
        }

        match self.calc_mean_and_var_of_matrix(res_group, bias.view()) {
            Some(((scale, _), _)) => {
                // println!("{}", scale);
                bias.map_inplace(|x| if *x == 0.0 {*x = f64::NAN} else { *x /= scale.sqrt()});
            },
            _ => {
                println!("Problem with computing mean. Skip scaling.");
                return None;
            }
        }
        Some(bias)
    }


    fn calc_mean_and_var_of_matrix(&self, res_group: &ResGroup, bias: ArrayView1<f64>) -> Option<((f64, f64), Array1<f64>)> {
        let mut res = Array1::<f64>::zeros((res_group.get_n_bins(),));
        for (bins1, bins2, counts) in res_group.get_raw_pixel_iter(CHUNKSIZE) {
            let data = self.pipe_product(res_group.get_n_bins(), bias, bins1, bins2, counts);
            res += &data;
        }

        let nnz_inds = res.mapv(|x| x != 0.0);
        if nnz_inds.is_empty() { return None; }
        let nnz_elems = Array1::from(utils::get_array_wrt_predicate(nnz_inds.view(), res.view()));
        nnz_elems.mean().zip(nnz_elems.central_moment(2).ok()).zip(Some(res))
    }


    fn pipe_binarize(&self, n_bins: usize, bins1: Array1<u32>, bins2: Array1<u32>, counts: Array1<u32>) -> Array1<u32> {
        let data = self.zeroing_diags(bins1.view(), bins2.view(), counts);
        let data = self.binarize(data);
        self.marginalize(n_bins, bins1.view(), bins2.view(), data.view())
    }

    fn pipe_zeroing(&self, n_bins: usize, bins1: Array1<u32>, bins2: Array1<u32>, counts: Array1<u32>) -> Array1<u32> {
        let data = self.zeroing_diags(bins1.view(), bins2.view(), counts);
        self.marginalize(n_bins, bins1.view(), bins2.view(), data.view())
    }

    fn pipe_product(&self, n_bins: usize, bias: ArrayView1<f64>, bins1: Array1<u32>, bins2: Array1<u32>, counts: Array1<u32>) -> Array1<f64> {
        let data = self.zeroing_diags(bins1.view(), bins2.view(), counts)
            .mapv(|x| x as f64);
        let data = self.outer_product(bias, bins1.view(), bins2.view(), data);
        self.marginalize(n_bins, bins1.view(), bins2.view(), data.view())
    }

    fn zeroing_diags(&self, bins1: ArrayView1<u32>, bins2: ArrayView1<u32>, mut counts: Array1<u32>) -> Array1<u32> {
        Zip::from(&mut counts).and(bins1).and(bins2).par_apply(|c, b1, b2| {
            let diff = if b1 > b2 {b1 - b2} else {b2 - b1};
            if diff < self.ignore_diags { *c = 0 };
        });
        counts
    }

    fn binarize(&self, mut data: Array1<u32>) -> Array1<u32> {
        data.map_mut(|x| if *x != 0 {*x = 1});
        data
    }

    fn marginalize<T>(&self, n_bins: usize, bins1: ArrayView1<u32>, bins2: ArrayView1<u32>, data: ArrayView1<T>)
        -> Array1<T> where T: Copy + ops::AddAssign + identities::Zero {
        let m1 = utils::bincount(n_bins, bins1, data);
        let m2 = utils::bincount(n_bins, bins2, data);
        m1 + m2
    }

    fn outer_product(&self, bias: ArrayView1<f64>, bins1: ArrayView1<u32>, bins2: ArrayView1<u32>, data: Array1<f64>) -> Array1<f64> {
        let bin1_bias = Array1::from_iter(bins1.iter().map(|&i| { bias[i as usize] }));
        let bin2_bias = Array1::from_iter(bins2.iter().map(|&i| { bias[i as usize] }));
        data * (bin1_bias * bin2_bias)
    }

    fn filter_by_predicate(&self, predicate: ArrayView1<bool>, mut array: Array1<f64>) -> Array1<f64> {
        azip!((array in &mut array, predicate in predicate) if *predicate {*array = 0.0});
        array
    }
}

fn debug(bias: &Array1<f64>, file: &Path) {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    let f = File::create(file).expect("Problem with file");
    let mut f = BufWriter::new(f);
    for x in bias.view() {
        writeln!(f, "{}", *x).expect("Problem with writing file");
    }
    f.flush().expect("Problem with flushing");
}