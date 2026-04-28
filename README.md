# A C++ Implementation for "Attacks on Goldreich's Pseudorandom Generators by Grouping and Solving"

This artifact contains source code for:

- full attack pipeline (`output_generator` + `attack_program`);
- toy end-to-end attack (`toy_attack`);
- FiLIP-oriented statistical verification (`FiLIP_verifier`).

## Dependencies

Required:

- C++17 compiler (`g++` tested with MinGW GCC 8.2.0)
- CMake >= 3.16 (optional but recommended)

Optional:

- OpenMP-capable toolchain (for faster `attack_program`)
- PowerShell (for harness scripts under `scripts/`)

`attack_program` is portable across Windows, Linux, and macOS:

- Windows uses memory-mapped file I/O in `verify.cpp`
- Linux/macOS use POSIX `mmap` in `verify.cpp`

## Platform Notes

The included scripts and commands were validated on:

- Windows + PowerShell
- MinGW GCC 8.2.0

The implementation is no longer Windows-only. It can also be built and run in standard Unix-like environments, including WSL, Linux, and macOS.

## Build

### Recommended (CMake)

```sh
cmake -S . -B build
cmake --build build
```

Generated binaries:

- `FiLIP_verifier`
- `toy_attack`
- `output_generator`
- `attack_program`

### Manual build

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

## Run Instructions

### 1) Statistical verification

```sh
./FiLIP_verifier
```

Expected output includes a table with `Good Group` and `Bad Group` probabilities.

How to interpret this output:

- `Good Group` and `Bad Group` are the main statistical quantities of interest.
- A successful sanity check is that the `Good Group` probability is clearly stronger than the `Bad Group` probability.
- For a faithful reproduction, these reported probabilities should also be reasonably close to the corresponding probability values reported in the paper's table for this parameter setting.
- `SAMPLE_SIZE` in `FiLIP_verify.cpp` controls the Monte Carlo accuracy of this estimate.
- Using a larger `SAMPLE_SIZE` (eg. 100000000) usually makes the observed probabilities more stable and closer to the paper-level values, at the cost of a longer runtime.
- This executable is the artifact-side statistical check corresponding to the FiLIP-style grouping behavior discussed in the paper.

### 2) Toy attack

```sh
./toy_attack
```

Expected output includes either:

- `[SUCCESS] Key Recovered ...`
- or `[FAIL] Could not recover key ...`

How to interpret this output:

- This is a small end-to-end demonstration of the attack logic on reduced parameters.
- It is intended as a fast correctness check for the attack workflow, not as a paper-scale benchmark.
- Because it is randomized, repeated runs may either succeed or fail.

### 3) Full attack

Step A: generate PRG outputs:

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

In practice, you should reserve at least `650 GB` of free disk space before running the default full dataset generation.

Runtime note:

- the full default generation can take many hours
- the exact time depends strongly on CPU speed, storage bandwidth, and whether the output directory is on SSD, HDD, a native Linux filesystem, or a mounted Windows filesystem
- if you want a machine-specific estimate first, reduce `num_instances` in `output_make.cpp`, run a smaller pilot job, and extrapolate from the observed throughput

Step B: launch attack:

```sh
./attack_program
```

Expected runtime output includes:

- selected `t_indices`
- number of filtered equations
- recovered-key verification status on fixed output bits

How to interpret this output:

- `selected t_indices` identifies the guessed group being tested
- `number of filtered equations` indicates how many generated samples survived the grouping filter for that guess
- `Recovered key passed PRG-output verification ...` is the main success criterion for the full artifact: it means the candidate key also matches independently generated PRG output bits, not only the linear system
- if no candidate passes the final verification step, then that attempted group did not yield a successful recovery within the configured number of trials

Relation to the paper:

- `FiLIP_verifier` corresponds to the statistical grouping behavior discussed in the paper
- `toy_attack` is only a reduced demonstration for quick validation
- `output_generator` + `attack_program` together correspond to the full grouped-equation recovery workflow described in the paper

## Configuration / Modes

The artifact supports two practical usage modes:

1. `quick-check` mode:
- run only `FiLIP_verifier` and `toy_attack`
- no large dataset generation required

2. `full-attack` mode:
- run `output_generator` then `attack_program`
- uses on-disk grouped equations and full recovery workflow

Key runtime constants can be changed directly in source files:

- `output_make.cpp`: dataset size (`num_instances`, `instances_per_file`)
- `attack.cpp`: number of inner solve attempts and verification bits
- `FiLIP_verify.cpp`: statistical experiment parameters

## Parameter Sets and Table 4

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

The relevant files are:

- `FiLIP_verify.cpp`
- `toy.cpp`
- `output_make.cpp`
- `attack.cpp`

## Test Harness and Summary Scripts

Run full harness:

```powershell
.\scripts\run_artifact_harness.ps1
```

This script:

- builds binaries
- runs `FiLIP_verifier`
- runs `toy_attack`
- writes logs to `artifact-logs/`
- generates summary via `scripts/summarize_harness_logs.ps1`

You can re-run summary only:

```powershell
.\scripts\summarize_harness_logs.ps1 -LogsDir artifact-logs
```

## Validation Status

The current codebase was re-checked after the Unix `mmap` update in `verify.cpp`.

Verified locally:

- all four binaries compile successfully with `g++ -std=c++17 -O3`
- `FiLIP_verifier` runs to completion and prints the expected probability table
- `toy_attack` runs in quick-check mode; because it is randomized, repeated runs may either succeed or fail
- `output_generator` and `attack_program` were validated at compile time, but the full end-to-end run was not repeated here because the default dataset generation writes about `650 GB`

## Source Code Organization

- `prng.h`, `prng.cpp`:
  PRG implementation and XOR-MAJ filter.
- `output_make.cpp`:
  dataset generation (`data/output_*.bin`).
- `verify.h`, `verify.cpp`:
  grouping filters, matrix construction, GF(2) solving, and output-based verification helpers.
- `attack.cpp`:
  full attack driver.
- `toy.cpp`:
  small parameter end-to-end demonstration.
- `FiLIP_verify.cpp`:
  high-fidelity statistical verification for FiLIP-style setting.
- `scripts/`:
  harness and summary scripts for artifact evaluation.

## License

This repository is released under the MIT License. See [LICENSE](LICENSE).
