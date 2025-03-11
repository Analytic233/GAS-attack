# Bit-fixing Correlation Attack on Random Local Functions

This is a simple implementation of the Bit-fixing Correlation Attack on Random Local Functions based on C++.

## Features and Future Improvements
Due to time constraints, the following optimization techniques have not yet been implemented:
- **Parallel Computing**: This will help in speeding up the computation process.
- **Sparse Matrix Gaussian Elimination**: To efficiently solve large, sparse systems.

These features will be added in the future.

## Key Information
The key used in this implementation is derived from hashing the value **"CRYPTO2025 "** using the standard **SHA-256** function.

## Usage Instructions
To run the program, you need to generate outputs first and then perform the attack. Here are the commands to do both steps in sequence:

```sh
# Step 1: Generate Outputs
g++ -std=c++17 prng.cpp output_make.cpp -o output_generator
./output_generator

# Step 2: Perform Attack
g++ -std=c++17 attack.cpp prng.cpp verify.cpp -o attack_program
./attack_program
