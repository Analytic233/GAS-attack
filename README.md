# A C++ Implementation for "Attacks on Goldreich's Pseudorandom Generators by Grouping and Solving"

This repository contains the C++ source code for the cryptanalysis presented in our paper. It includes a full-scale attack implementation, a self-contained toy example, and dedicated verification scripts designed to validate the core statistical claims of our method, particularly in the context of FiLIP-like constructions.

## Getting Started

We recommend starting with the verification scripts, as they are fast and directly demonstrate the core principles of our attack.

### 1. Quick Verification & Demonstration (Recommended for Reviewers)

These scripts are designed to run quickly on a standard laptop and provide immediate insight into the attack's mechanics.

#### A. Statistical Verification Script (`Filp_verify.cpp`)

This script validates the theoretical probabilities presented in our paper. It uses a "fakely done by hand" approach, to analyze the statistical behavior of different groups without performing a full search.

**Compile and Run:**
```sh
# Compile with optimizations for fast execution
g++ -std=c++17 -O3 Filp_verify.cpp -o filp_verifier

# Run the verification
./filp_verifier
```

#### B. Self-Contained Toy Example (`toy.cpp`)

This is a complete, end-to-end attack on a scaled-down parameter set. It demonstrates the entire attack pipeline, from grouping to final key recovery.

**Compile and Run:**
```sh
# Compile the toy example
g++ -std=c++17 -O3 toy.cpp -o toy_attack

# Run the attack
./toy_attack
```
This program will recover the secret key for the toy instance, typically completing in less than a minute.

### 2. Full-Scale Attack on the PRG Instance ($N=256$)

This is the implementation of the primary attack presented in our paper. **Warning:** This process is computationally intensive and requires significant time and disk space.

#### Step 1: Generate PRG Output Samples

First, generate a large dataset of PRG outputs.

```sh
# Compile the data generator
g++ -std=c++17 -O3 prng.cpp output_make.cpp -o output_generator

# Run the generator (this may take several hours and generate gigabytes of data)
./output_generator
```

#### Step 2: Perform the Attack

With the data generated, run the main attack program to recover the secret seed.

```sh
# Compile the main attack program with optimizations and OpenMP for parallel acceleration
g++ -std=c++17 -O3 -fopenmp attack.cpp prng.cpp verify.cpp -o attack_program

# Run the attack (this is the most time-consuming step)
./attack_program
```

## Key Derivation

For reproducibility, the default secret key used in all implementations is the **SHA-256** hash of the string **"EUROCRYPTO2026"**.

