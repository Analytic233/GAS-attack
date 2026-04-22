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

`attack_program` is now portable across Windows, Linux, and macOS:

- Windows uses memory-mapped file I/O in `verify.cpp`
- Linux/macOS use POSIX `mmap` in `verify.cpp`

## Platform Used For Reported Runs

The included scripts and commands were validated on:

- Windows + PowerShell
- MinGW GCC 8.2.0

## Build

### Recommended (CMake)

```sh
cmake -G "MinGW Makefiles" -S . -B build
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

### 2) Toy attack

```sh
./toy_attack
```

Expected output includes either:

- `[SUCCESS] Key Recovered ...`
- or `[FAIL] Could not recover key ...`

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

Step B: launch attack:

```sh
./attack_program
```

Expected runtime output includes:

- selected `t_indices`
- number of filtered equations
- recovered-key verification status on fixed output bits

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

## Test Harness and Summary Scripts

Run full harness:

```powershell
.\scripts\run_artifact_harness.ps1
```

This script:

- builds binaries;
- runs `FiLIP_verifier`;
- runs `toy_attack`;
- writes logs to `artifact-logs/`;
- generates summary via `scripts/summarize_harness_logs.ps1`.

You can re-run summary only:

```powershell
.\scripts\summarize_harness_logs.ps1 -LogsDir artifact-logs
```

## Validation Status

The current codebase was re-checked after the Linux `mmap` update in `verify.cpp`.

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

## Interpreting Main Outputs

- `FiLIP_verifier`:
  inspect final `Good Group` / `Bad Group` probability lines.
- `toy_attack`:
  check for success/fail line and recovered key print.
- `attack_program`:
  success means candidate key passed fixed-bit output comparison; failure means no candidate passed within configured trials.

## License

This repository is released under the MIT License. See [LICENSE](LICENSE).
