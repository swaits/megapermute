# megapermute

A million permutation hypothesis test utility.

## Installation

1. [install Rust](https://www.rust-lang.org/tools/install)
2. `cargo install --git https://git.sr.ht/~swaits/megapermute`

## Usage

`megapermute` requires two input files, `control.dat` and `treatment.dat`. Each is a text
file containing a list of numbers, treated internally as `f64`.

- `control.dat` the control sample
- `treatment.dat` the treatment sample

With these two files in the current directly, running `megapermute` computes
the empircal difference in means and p-value for the null hypothesis (H_0 ==
H_a).
