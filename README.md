# Eurocrypt Artifact for "Attacks on Goldreich's Pseudorandom Generators by Grouping and Solving"

This is the artifact belonging to the paper:

> Ximing Fu, Mo Li, Shihan Lyu, and Chuanyi Liu, Attacks on Goldreich's Pseudorandom Generators by Grouping and Solving, Eurocrypt 2026.

The repository contains:

1. C++ implementations for the main attack workflow, including `output_generator`, `attack_program`, `toy_attack`, and `FiLIP_verifier`.
2. Utility scripts in `scripts/` for artifact harnessing and for computing Table 5 and Table 6 values directly from the formulas in the paper.

# Dependencies

Required dependencies:

- C++17 compiler (`g++` tested with MinGW GCC 8.2.0)
- CMake >= 3.16

Optional dependencies:

- OpenMP-capable toolchain (for faster `attack_program`)
- PowerShell (for the harness scripts in `scripts/`)
- Python 3 (for the Table 5 / Table 6 calculation scripts in `scripts/`)

The attack code is portable across Windows, Linux, and macOS:

- Windows uses memory-mapped file I/O in `verify.cpp`
- Linux/macOS use POSIX `mmap` in `verify.cpp`

The included commands were validated locally on:

- Windows + PowerShell
- MinGW GCC 8.2.0

The implementation can also be built and run in standard Unix-like environments, including WSL, Linux, and macOS.

# Running the Artifact

The artifact supports two complementary workflows:

1. Running the C++ programs for statistical verification, the toy attack, and the full grouped-equation attack pipeline.
2. Running Python scripts that compute Table 5 and Table 6 values directly from the paper's formulas.

# Build

## Recommended (CMake)

```sh
cmake -S . -B build
cmake --build build
```

This builds:

- `FiLIP_verifier`
- `toy_attack`
- `output_generator`
- `attack_program`

## Manual build

```sh
g++ -std=c++17 -O3 FiLIP_verify.cpp -o FiLIP_verifier
g++ -std=c++17 -O3 toy.cpp -o toy_attack
g++ -std=c++17 -O3 prng.cpp output_make.cpp -o output_generator
g++ -std=c++17 -O3 attack.cpp prng.cpp verify.cpp -o attack_program
```

If your toolchain supports OpenMP, you may use:

```sh
g++ -std=c++17 -O3 -fopenmp attack.cpp prng.cpp verify.cpp -o attack_program
```

# Quick Check

For a fast sanity check of the artifact, run:

```sh
./FiLIP_verifier
./toy_attack
```

Interpretation:

- `FiLIP_verifier` is the statistical validation executable for the FiLIP-style grouping behavior discussed in the paper.
- `toy_attack` is a reduced end-to-end demonstration of the attack logic on small parameters.
- `toy_attack` is randomized, so repeated runs may either succeed or fail.

# Full Attack Workflow

## Step A: Generate PRG outputs

```sh
./output_generator
```

This creates `data/output_*.bin` files.

Default storage usage of `output_generator`:

- each sample stores `1 + 32 + 32 = 65` bytes
- default `num_instances = 10,000,000,000`
- default `instances_per_file = 20,000,000`
- number of output files is `500` (`10,000,000,000 / 20,000,000`)
- each output file is `1,300,000,000` bytes, about `1.30 GB` or `1.21 GiB`
- total output size is `650,000,000,000` bytes, about `650 GB` or `605.36 GiB`

In practice, reserve at least `650 GB` of free disk space before running the default full dataset generation.

Runtime note:

- the full default generation can take many hours
- the exact runtime depends strongly on CPU speed and storage bandwidth
- for a machine-specific estimate, reduce `num_instances` in `output_make.cpp`, run a smaller pilot job, and extrapolate from the observed throughput

## Step B: Launch attack

```sh
./attack_program
```

Expected output includes:

- selected `t_indices`
- number of filtered equations
- recovered-key verification status on fixed output bits

Interpretation:

- `selected t_indices` identifies the guessed group being tested
- `number of filtered equations` indicates how many generated samples survived the grouping filter for that guess
- `Recovered key passed PRG-output verification ...` is the main success criterion for the full artifact
- if no candidate passes the final verification step, then the attempted group did not yield a successful recovery within the configured number of trials

Relation to the paper:

- `FiLIP_verifier` corresponds to the statistical grouping behavior discussed in the paper
- `toy_attack` is only a reduced demonstration for quick validation
- `output_generator` + `attack_program` together correspond to the full grouped-equation recovery workflow described in the paper

# Statistical Verification

The statistical verification executable is:

```sh
./FiLIP_verifier
```

Expected output includes a table with `Good Group` and `Bad Group` probabilities.

Interpretation:

- `Good Group` and `Bad Group` are the main statistical quantities of interest
- a successful sanity check is that the `Good Group` probability is clearly stronger than the `Bad Group` probability
- for a faithful reproduction, these reported probabilities should also be reasonably close to the corresponding probability values reported in the paper's table for this parameter setting
- `SAMPLE_SIZE` in `FiLIP_verify.cpp` controls the Monte Carlo accuracy of this estimate
- using a larger `SAMPLE_SIZE` usually makes the observed probabilities more stable and closer to the paper-level values, at the cost of a longer runtime

# Parameter Sets and Table 4

The current codebase hardcodes one representative parameter set for the included artifact runs. It does **not** currently expose all parameter sets from Table 4 of the paper as command-line options or a unified configuration layer.

For readers comparing the code to Table 4:

- the artifact demonstrates the workflow on one concrete parameter set
- additional Table 4 parameter sets are not yet packaged as separate ready-to-run presets in the current repository
- reproducing other Table 4 rows currently requires manual editing of the constants in the relevant source files

Table 4 seeds and their `SHA-256` keys:

- `SHA256(EUROCRYPTO2025)`:
  `8b0bac2acc717bb466eaabed78e2fbdb3ad98960debdc5532b67ec5aa665546f`
- `SHA256(CRYPTO2025)`:
  `85e679140bf62c020f22574f8c887cef4eb86c24c962f3e7d12e22a28c078fb1`
- `SHA256(EUROCRYPT2026)`:
  `b72412e79157cf37cbfbe7ee7b5b02d9cf5096253efa99e4863e78d8055cd6f5`

The currently hardcoded artifact key in `attack.cpp` / `output_make.cpp` corresponds to:

- `SHA256(EUROCRYPT2026)`:
  `0xb7, 0x24, 0x12, 0xe7, 0x91, 0x57, 0xcf, 0x37, 0xcb, 0xfb, 0xe7, 0xee, 0x7b, 0x5b, 0x02, 0xd9, 0xcf, 0x50, 0x96, 0x25, 0x3e, 0xfa, 0x99, 0xe4, 0x86, 0x3e, 0x78, 0xd8, 0x05, 0x5c, 0xd6, 0xf5`

# Formula-Based Table Scripts

The repository also includes scripts that compute paper-style table outputs for the FiLIP analysis without running the full experimental attack. These scripts are meant to help readers inspect the quantities reported in the paper through direct computation rather than through a long data-generation experiment.

## Table 5

```sh
python scripts/compute_table5_filip_xor_thr.py
```

This prints a paper-style table for FiLIP instantiated with `XORa-THRt,b`.

Output columns:

- `p`
- `Teg`
- `Ts`
- `Tbit`
- `Tenc`

These values are computed directly from the formulas used in the paper's analysis .


## Table 6

```sh
python scripts/compute_table6_filip_xor_thr_thr.py
```

This prints a paper-style table for FiLIP instantiated with `XORa-THRt1,b1-THRt2,b2`.

Output columns:

- `p`
- `D`
- `e`
- `ep`
- `Teg`
- `Ts`
- `Tbit`
- `Tenc`


# Harness Scripts

Run the local artifact harness with:

```powershell
.\scripts\run_artifact_harness.ps1
```

This script:

- builds the binaries
- runs `FiLIP_verifier`
- runs `toy_attack`
- writes logs to `artifact-logs/`
- generates a short summary via `scripts/summarize_harness_logs.ps1`

You can re-run the summary only with:

```powershell
.\scripts\summarize_harness_logs.ps1 -LogsDir artifact-logs
```

# Validation Status

Verified locally:

- all four main binaries compile successfully with `g++ -std=c++17 -O3`
- `FiLIP_verifier` runs to completion and prints the expected probability table
- `toy_attack` runs in quick-check mode, though repeated runs may either succeed or fail because it is randomized
- `output_generator` and `attack_program` were validated at compile time, but the full default end-to-end run was not repeated here because the default dataset generation writes about `650 GB`

# Organization of Files

```text
.
|-- attack.cpp
|-- CMakeLists.txt
|-- FiLIP_verify.cpp
|-- LICENSE
|-- output_make.cpp
|-- prng.cpp
|-- prng.h
|-- README.md
|-- toy.cpp
|-- verify.cpp
|-- verify.h
`-- scripts
    |-- compute_table5_filip_xor_thr.py
    |-- compute_table6_filip_xor_thr_thr.py
    |-- run_artifact_harness.ps1
    `-- summarize_harness_logs.ps1
```

# License

This repository is released under the MIT License. See [LICENSE](LICENSE).
