use std::fs::File;
use std::path::Path;
use serde::Deserialize;
use ascii::{AsciiString, AsAsciiStr};
use std::error::Error;
use std::cmp::Ordering;
use ndarray::{azip, Array1, ArrayView1, self};
use num_traits::identities;
use std::ops;

// pub const CHUNKSIZE: usize = 50_000_000;

#[derive(Debug, Deserialize)]
struct Record<'a> {
    tig_name: &'a str,
    length: u64,
}

pub fn parse_tig_lengths(file_name: &Path) -> Result<Vec<(AsciiString, u64)>, Box<dyn Error>> {
    let mut tig_lengths: Vec<(AsciiString, u64)> = Vec::new();
    let file = File::open(file_name)?;

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(file);
    let mut raw_record = csv::ByteRecord::new();

    while rdr.read_byte_record(&mut raw_record)? {
        let record: Record = raw_record.deserialize(None)?;
        let nm = AsciiString::from(record.tig_name.as_ascii_str()?);
        tig_lengths.push((nm, record.length));
    }

    tig_lengths.sort_by_key(|x| x.1);
    Ok(tig_lengths)
}

// https://rosettacode.org/wiki/Quickselect_algorithm#Rust
pub fn get_array_wrt_predicate<T: Copy>(predicate: ArrayView1<bool>, array: ArrayView1<T>) -> Vec<T> {
    assert_eq!(predicate.len(), array.len());
    ndarray::Zip::from(predicate).and(array).fold(Vec::<T>::new(), |mut v, &p, &a| {
        if p { v.push(a) };
        v
    })
}

pub fn get_vec_wrt_predicate<T: Copy>(predicate: &[bool], array: &[T]) -> Vec<T> {
    assert_eq!(predicate.len(), array.len());
    predicate.iter().zip(array.iter()).fold(Vec::<T>::new(), |mut v, (&p, &a)| {
        if p { v.push(a) }; v
    })
}


pub fn get_zooming_order(resolutions: &[u32]) -> Vec<i32> {
    let mut pred = vec![-1; resolutions.len()];

    for (i, res) in resolutions.iter().enumerate().skip(1) {
        let mut p = (i - 1) as i32;

        while p >= 0 {
            if res % resolutions[p as usize] == 0 {
                pred[i] = p;
                break;
            }
            p -= 1;
        }
    }

    pred
}

pub fn bincount<Q>(length: usize, array: ArrayView1<u32>, weights: ArrayView1<Q>) -> Array1<Q>
    where Q: Copy + ops::AddAssign + identities::Zero {
    let mut counts: Array1<Q> = Array1::zeros((length,));

    if array.len() != weights.len() {
        eprintln!("Houston, we have a problem.");
        //TODO
        return counts;
    }

    azip!((&array in array, &weights in weights) counts[array as usize] += weights);
    counts
}

pub fn mad(mut data: Vec<f64>) -> Option<f64>{
    median(&data).and_then(|med| {
        data.iter_mut().for_each(|x| { *x = (*x - med).abs(); });
        median(&data)
    })
}

pub fn median(data: &[f64]) -> Option<f64> {
    data.iter().for_each(|x| {
        if !x.is_finite() {panic!("Try to find median of vector with nans or infinities")}
    });

    let size = data.len();
    match size {
        0 => None,
        even if even % 2 == 0 => {
            let fst_med = select(data, (even / 2) - 1);
            let snd_med = select(data, even / 2);

            match (fst_med, snd_med) {
                (Some(fst), Some(snd)) => Some((fst + snd) as f64 / 2.0),
                _ => None
            }
        },
        odd => select(data, odd / 2).map(|x| x as f64)
    }
}

fn select(data: &[f64], k: usize) -> Option<f64> {
    let part = partition(data);

    match part {
        None => None,
        Some((left, pivot, right)) => {
            let pivot_idx = left.len();

            match pivot_idx.cmp(&k) {
                Ordering::Equal => Some(pivot),
                Ordering::Greater => select(&left, k),
                Ordering::Less => select(&right, k - (pivot_idx + 1)),
            }
        },
    }
}


fn partition(data: &[f64]) -> Option<(Vec<f64>, f64, Vec<f64>)> {
    match data.len() {
        0 => None,
        _ => {
            let (pivot_slice, tail) = data.split_at(1);
            let pivot = pivot_slice[0];
            let (left, right) = tail.iter()
                .fold((vec![], vec![]), |mut splits, next| {
                    {
                        let (ref mut left, ref mut right) = &mut splits;
                        if next < &pivot {
                            left.push(*next);
                        } else {
                            right.push(*next);
                        }
                    }
                    splits
                });

            Some((left, pivot, right))
        }
    }
}
