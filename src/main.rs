use std::{
    fs::File,
    io::{BufRead, BufReader},
    iter,
};

use anyhow::Result;
use rand::{prelude::SliceRandom, thread_rng};
use rayon::prelude::*;

// the number of threads to pass to rayon (which smartly splits it into actual number of threads)
const N_THREADS: usize = 1_000;

// the number of permutation tests to run per "thread" (above)
const N_PERMUTATIONS_PER_THREAD: usize = 1_000;

fn load_f64s(filename: &str) -> Result<Vec<f64>> {
    let f = BufReader::new(File::open(filename)?);
    Ok(f.lines()
        .map(|l| {
            l.expect("input error: unable to read line")
                .parse::<f64>()
                .expect("parse error: unable to parse f64")
        })
        .collect())
}

// accepts an iterator of f64's and computes mean using Welford's online algorithm
fn mean<'a>(iter: impl Iterator<Item = &'a f64>) -> f64 {
    iter.enumerate()
        .fold(0.0, |mu, (i, x)| mu + ((x - mu) / (i + 1) as f64))
}

#[derive(Clone, PartialEq, Eq)]
enum Group {
    Control,
    Treatment,
}

fn permute(control: &[f64], treatment: &[f64], mu_diff: f64) -> f64 {
    let count: f64 = (0..N_THREADS)
        .into_par_iter()
        .map(|_| {
            // create an rng so we can shuffle later
            let mut rng = thread_rng();

            // create an index array of [Control, Treatment] which we'll shuffle repeatedly to make our
            // selections during each permutation
            let mut index: Vec<Group> = iter::repeat(Group::Control)
                .take(control.iter().len())
                .chain(iter::repeat(Group::Treatment).take(treatment.iter().len()))
                .collect();

            // do the actual permutations
            (0..N_PERMUTATIONS_PER_THREAD)
                .filter(|_| {
                    // shuffle our group selections
                    index.shuffle(&mut rng);

                    // comput emean of this "control permutation"
                    let permuted_mean_control = mean(
                        control
                            .iter()
                            .chain(treatment.iter())
                            .enumerate()
                            .filter_map(|(i, x)| {
                                if index[i] == Group::Control {
                                    Some(x)
                                } else {
                                    None
                                }
                            }),
                    );

                    // comput emean of this "treatment permutation"
                    let permuted_mean_treatment = mean(
                        control
                            .iter()
                            .chain(treatment.iter())
                            .enumerate()
                            .filter_map(|(i, x)| {
                                if index[i] == Group::Treatment {
                                    Some(x)
                                } else {
                                    None
                                }
                            }),
                    );

                    // select this permutation if the diff in these permuted means is more than the
                    // empirical diff of means
                    (permuted_mean_treatment - permuted_mean_control) > mu_diff
                })
                .count() as f64
        })
        .sum();

    // p-value is the ratio of permutations where delta(mean) exceeded empircal delta(mean) to
    // total permutations
    let p_value = count / (N_THREADS * N_PERMUTATIONS_PER_THREAD) as f64;

    // adjust for left or right tail
    if mu_diff < 0.0 {
        1.0 - p_value
    } else {
        p_value
    }
}

fn main() -> Result<()> {
    // read data
    let control = load_f64s("control.dat")?;
    let treatment = load_f64s("treatment.dat")?;

    // compute empircal means
    let mean_control = mean(control.iter());
    let mean_treatment = mean(treatment.iter());
    println!("                 mu_control = {}", mean_control);
    println!("               mu_treatment = {}", mean_treatment);
    println!(
        "(mu_treatment - mu_control) = {}",
        mean_treatment - mean_control
    );

    // run permutation test to compute p-value
    let pvalue = permute(&control, &treatment, mean_treatment - mean_control);
    println!("                    p-value = {}", pvalue);

    Ok(())
}

#[cfg(test)]
mod tests {
    // import names from main scope
    use super::*;

    const MEAN_EPSILON: f64 = 0.000001;
    const PVALUE_EPSILON: f64 = 0.001;

    // for float comparison
    fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    #[test]
    fn test_mouse_data() {
        // test data from Table 2.1 in "An Introduction to the Bootstrap" (book)
        let control = vec![52.0, 104.0, 146.0, 10.0, 51.0, 30.0, 40.0, 27.0, 46.0];
        let treatment = vec![94.0, 197.0, 16.0, 38.0, 99.0, 141.0, 23.0];

        // compute empircal means
        let mean_control = mean(control.iter());
        assert!(approx_eq(mean_control, 56.22222222222222, MEAN_EPSILON));
        let mean_treatment = mean(treatment.iter());
        assert!(approx_eq(mean_treatment, 86.85714285714286, MEAN_EPSILON));
        assert!(approx_eq(
            mean_treatment - mean_control,
            30.63492063492064,
            MEAN_EPSILON
        ));

        // run permutation test to compute p-value
        let pvalue = permute(&control, &treatment, mean_treatment - mean_control);
        assert!(approx_eq(pvalue, 0.13896357, PVALUE_EPSILON));
    }
}
