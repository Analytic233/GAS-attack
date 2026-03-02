#ifndef VERIFIER_UTILS_H
#define VERIFIER_UTILS_H

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#define MAX_SIZE 256

// Convert key bytes to a bit string (MSB first per byte).
std::string bytesToBitstring(const std::vector<uint8_t> &byte_data);
// Convert a '0'/'1' bit string to bytes.
std::vector<uint8_t> bitstringToBytes(const std::string &bitstring);

// Build deterministic challenge subsets for reproducible output-based verification.
std::vector<std::vector<int>> generateDeterministicSubsets(int num_samples, int N = 256, int n = 74, uint32_t seed = 20260228u);
// Evaluate PRG outputs for one key on the provided challenge subsets.
std::vector<uint8_t> generateOutputsFromKey(const std::vector<uint8_t> &key_bytes, const std::vector<std::vector<int>> &subsets, int a = 10, int b = 64);
// Verify candidate key by comparing outputs against expected reference outputs.
bool verifyCandidateByPRGOutputs(const std::string &candidate_key_bits, const std::vector<std::vector<int>> &subsets, const std::vector<uint8_t> &expected_outputs, int a = 10, int b = 64);

// XOR-MAJ filter function.
int xorMaj(const std::vector<int> &xor_bits, const std::vector<int> &maj_bits);

// Generate one subset/permutation that includes t_indices in the MAJ segment.
std::vector<int> generatePermutedIndices(const std::vector<int> &t_indices, int total_size = 74, int total_range = 256);

// Recover full key bit-string from solved variables + fixed t_indices.
std::string recover_key(const std::vector<int> &solution, const std::vector<int> &t_indices, int total_key_size = 256);

void printMatrix(const std::vector<std::vector<int>> &matrix);

// Build A*x=b from grouped instances (at most 245 rows sampled each attempt).
std::pair<std::vector<std::bitset<MAX_SIZE>>, std::vector<int>> generateMatrixFromInstances(
    const std::vector<std::pair<int, std::vector<int>>> &successfulInstances, const std::vector<int> &tIndices, int totalKeySize = 256);

// Keep only solutions satisfying A*x=b.
std::vector<std::vector<int>> verifySolutions(const std::vector<std::bitset<MAX_SIZE>> &A, const std::vector<int> &b, const std::vector<std::vector<int>> &solutions);

bool checkMajPart(const std::vector<int> &guess_indices, const std::vector<int> &permuted_indices);

// Solve linear system over GF(2) with optional free-variable expansion.
std::vector<std::vector<int>> gaussianEliminationMod2WithFreeVariables(std::vector<std::bitset<MAX_SIZE>> &A, std::vector<int> &b, int rows, int cols);

std::vector<int> selectRandomIndices(int t);

// Hamming weight of XOR between two byte arrays.
int hammingWeight(const uint8_t *a, const uint8_t *b, size_t size = 32);

void writeResultToFile(std::ofstream &outFile, uint8_t result_bit, const uint8_t *xor_vector, const uint8_t *maj_vector, size_t size = 32);
void processSingleFileWithMemoryMapping(const std::string &inputFilePath, const std::vector<int> &t_indices, std::ofstream &outFile);
void processFilesWithMemoryMapping(const std::string &folderPath, const std::vector<int> &t_indices, const std::string &outputFilePath, int num_files);

#endif // VERIFIER_UTILS_H
