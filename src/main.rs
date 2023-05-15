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

// load a file of numbers as f64 and panic if we run into any problems.
// This function takes a filename as a string and returns a vector of f64s.
// If there is an error, it will return a Result with an error message.
fn load_f64s(filename: &str) -> Result<Vec<f64>> {
    // Open the file and read it into a buffer.
    let f = BufReader::new(File::open(filename)?);
    // Map each line of the file to a f64.
    Ok(f.lines()
        .map(|l| {
            // If there is an error reading the line, panic.
            l.expect("input error: unable to read line")
                // If there is an error parsing the line, panic.
                .parse::<f64>()
                .expect("parse error: unable to parse f64")
        })
        // Collect the f64s into a vector.
        .collect())
}

// This function accepts an iterator of f64's and computes mean using Welford's online algorithm
fn mean<'a>(iter: impl Iterator<Item = &'a f64>) -> f64 {
    // The function uses the enumerate method to get the index and value of each element in the iterator
    // The fold method is used to iterate through the iterator and add the values to the accumulator
    // The accumulator is initialized to 0.0
    // The accumulator is updated by adding the difference between the current value and the accumulator divided by the index plus 1
    iter.enumerate()
        .fold(0.0, |mu, (i, x)| mu + ((x - mu) / (i + 1) as f64))
}

// this enum is used to denote group memebership during each permutation
#[derive(Clone, PartialEq, Eq)]
enum Group {
    Control,
    Treatment,
}

// Run the permutation tests.
//
// Accepts the control and treatment arrays along with the difference in empircal means (treatment
// minus control).
//
// Runs N_PERMUTATIONS_PER_THREAD on each of N_THREADS using Rayon's `par_iter()`.
//
// The original data (arrays) is never copied. Each of the iterations creates an index array of
// `enum Group` denoting `Control` or `Treatment` at each index. That `index` array is shuffled to
// permute group memebership. Finally, each of the group's means are computed and compared, and
// filtered to only count differences which are larger than the `mu_diff` parameter.
//
// The final p-value is the number of permutations where the diff in means exceeded `mu_diff`
// divided by (N_THREADS * N_PERMUTATIONS_PER_THREAD).
//
// Note: the left tail or right tail is chosen automatically based on whether the empircal mean
// delta is positive or negative.
fn permutation_test(control: &[f64], treatment: &[f64], mu_diff: f64) -> f64 {
    // use Rayon to divide this work across N_THREADS, ultimately counting the number of
    // permutations where delta(means) > mu_diff
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

                    // variables for mean computation
                    let (mut mu_control, mut n_control) = (0.0, 0.0);
                    let (mut mu_treatment, mut n_treatment) = (0.0, 0.0);

                    // walk the combined iterator and add each element to the corresponding mean
                    // using Welford's online algorithm
                    control.iter().chain(treatment.iter()).enumerate().for_each(
                        |(i, x)| match index[i] {
                            Group::Control => {
                                n_control += 1.0;
                                mu_control += (x - mu_control) / n_control
                            }
                            Group::Treatment => {
                                n_treatment += 1.0;
                                mu_treatment += (x - mu_treatment) / n_treatment
                            }
                        },
                    );

                    // select this permutation if the diff in these permuted means is more than the
                    // empirical diff of means
                    (mu_treatment - mu_control) > mu_diff
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

// convert p-value to conventional language
fn pvalue_to_string(p: f64) -> String {
    match p {
        p if p < 0.01 => "very strong evidence against null hypothesis",
        p if p < 0.025 => "strong evidence against null hypothesis",
        p if p < 0.05 => "reasonably strong evidence against null hypothesis",
        p if p < 0.10 => "borderline evidence against null hypothesis",
        _ => "no evidence against null hypothesis",
    }
    .to_string()
}

fn main() -> Result<()> {
    // read data
    let control = load_f64s("control.dat")?;
    let treatment = load_f64s("treatment.dat")?;

    // compute empircal means
    let mean_control = mean(control.iter());
    let mean_treatment = mean(treatment.iter());
    println!("                 mu_control = {}", mean_control);
    println!("                  N_control = {}", control.len());
    println!("               mu_treatment = {}", mean_treatment);
    println!("                N_treatment = {}", treatment.len());
    println!(
        "(mu_treatment - mu_control) = {}",
        mean_treatment - mean_control
    );

    // run permutation test to compute p-value
    let pvalue = permutation_test(&control, &treatment, mean_treatment - mean_control);
    println!("                    p-value = {}", pvalue);
    println!("                     result = {}", pvalue_to_string(pvalue));

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

    // This is a test file for the permutation test.
    // It tests the permutation test on the mouse data from Table 2.1 in "An Introduction to the Bootstrap" (book)
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
        let pvalue = permutation_test(&control, &treatment, mean_treatment - mean_control);
        assert!(approx_eq(pvalue, 0.13896357, PVALUE_EPSILON));
    }
}
