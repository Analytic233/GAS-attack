# A C++ Implementation for "Attacks on Goldreich's Pseudorandom Generators by Grouping and Solving"

This repository contains C++ code for the cryptanalysis in our paper, including:

- a full-scale attack implementation,
- a self-contained toy attack example,
- and verification scripts for key statistical claims (including FiLIP-like constructions).

## Build

### Recommended: CMake (Windows/Linux/macOS)

```sh
cmake -S . -B build
cmake --build build --config Release
```

If you are using MinGW on Windows, prefer:

```sh
cmake -G "MinGW Makefiles" -S . -B build
cmake --build build
```

Generated binaries:

- `FiLIP_verifier`
- `toy_attack`
- `output_generator`
- `attack_program`

### Manual g++ build (optional)

```sh
g++ -std=c++17 -O3 FiLIP_verify.cpp -o FiLIP_verifier
g++ -std=c++17 -O3 toy.cpp -o toy_attack
g++ -std=c++17 -O3 prng.cpp output_make.cpp -o output_generator
g++ -std=c++17 -O3 -fopenmp attack.cpp prng.cpp verify.cpp -o attack_program
# If your toolchain does not support OpenMP, remove -fopenmp.
```

## Quick Start

We recommend starting from fast verification scripts.

### 1. Statistical Verification (`FiLIP_verify.cpp`)

```sh
./FiLIP_verifier
```

### 2. End-to-End Toy Attack (`toy.cpp`)

```sh
./toy_attack
```

## Full-Scale Attack on PRG Instance (`N=256`)

This is computationally intensive and may require large disk space.

### Step 1: Generate PRG output samples

```sh
./output_generator
```

### Step 2: Run the main attack

```sh
./attack_program
```

## Key Derivation

For reproducibility, the default secret key is the SHA-256 hash of:

`EUROCRYPTO2026`
