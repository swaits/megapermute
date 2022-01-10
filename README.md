# megapermute

A million permutation hypothesis test utility.

## Installation

1. If you don't already have rust installed, [install it](https://www.rust-lang.org/tools/install).
2. Then run: `cargo install --git https://git.sr.ht/~swaits/megapermute`

## Usage

`megapermute` is a command-line utility which reads two input files,
`control.dat` and `treatment.dat`. Each is a text file containing a list of
numbers, treated internally as 64-bit floating point numbers.

- `control.dat`: the control sample
- `treatment.dat`: the treatment sample

With these two files in the current directory, running `megapermute` computes
the empircal difference in means and p-value for the null hypothesis, that H_0
is equivalent to H_a. You can think of p-value like the probability that the
null hypothesis is true.

### Example

Given the following two input files:

```bash
~/tmp
❯ bat control.dat
───────┬───────────────────────
       │ File: control.dat
───────┼───────────────────────
   1   │ 52
   2   │ 104
   3   │ 146
   4   │ 10
   5   │ 51
   6   │ 30
   7   │ 40
   8   │ 27
   9   │ 46
───────┴───────────────────────

~/tmp
❯ bat treatment.dat
───────┬───────────────────────
       │ File: treatment.dat
───────┼───────────────────────
   1   │ 94
   2   │ 197
   3   │ 16
   4   │ 38
   5   │ 99
   6   │ 141
   7   │ 23
───────┴──────────────────────
```

Running `megapermute` should yield results like this:

```bash
~/tmp
❯ megapermute
                 mu_control = 56.22222222222222
               mu_treatment = 86.85714285714286
(mu_treatment - mu_control) = 30.63492063492064
                    p-value = 0.139556
```

## Performance

My goal with this was to take advantage of multiple cores by using
[`rayon`](https://crates.io/crates/rayon) iterators to compute the permutation
hypothesis test quickly. As of January 2022, this application can do one
million permutations on the example dataset listed above in 39.4 ms (N=69,
SD=2.6 ms), which is on the order of 25M permutations/second.

```bash
~/tmp
❯ hyperfine megapermute
Benchmark 1: megapermute
  Time (mean ± σ):      39.4 ms ±   2.6 ms    [User: 499.3 ms, System: 10.2 ms]
  Range (min … max):    35.1 ms …  48.4 ms    69 runs
```

Benchmark machine info:

```bash
~/tmp
❯ neofetch
                    'c.          stephenwaits@SWaits-Mac-01.local
                 ,xNMM.          --------------------------------
               .OMMMMo           OS: macOS 11.6.1 20G224 x86_64
               OMMM0,            Host: MacBookPro16,1
     .;loddo:' loolloddol;.      Kernel: 20.6.0
   cKMMMMMMMMMMNWMMMMMMMMMM0:    Uptime: 31 days, 21 hours, 23 mins
 .KMMMMMMMMMMMMMMMMMMMMMMMWd.    Packages: 149 (brew)
 XMMMMMMMMMMMMMMMMMMMMMMMX.      Shell: zsh 5.8
;MMMMMMMMMMMMMMMMMMMMMMMM:       Resolution: 3840x2160, 1920x1080
:MMMMMMMMMMMMMMMMMMMMMMMM:       DE: Aqua
.MMMMMMMMMMMMMMMMMMMMMMMMX.      WM: Amethyst
 kMMMMMMMMMMMMMMMMMMMMMMMMWd.    Terminal: WezTerm
 .XMMMMMMMMMMMMMMMMMMMMMMMMMMk   CPU: Intel i9-9880H (16) @ 2.30GHz
  .XMMMMMMMMMMMMMMMMMMMMMMMMK.   GPU: Intel UHD Graphics 630, AMD Radeon Pro 5500M
    kMMMMMMMMMMMMMMMMMMMMMMd     Memory: 21514MiB / 32768MiB
     ;KMMMMMMMWXXWMMMMMMMk.
       .cooc,.    .,coo:.

```
