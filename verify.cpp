#include "verify.h"
#include "prng.h"
#include <algorithm>
#include <bitset>
#include <chrono>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <ctime>
#if defined(_MSC_VER)
#include <intrin.h>
#endif
#include <windows.h>

std::mt19937 generator(static_cast<unsigned int>(std::time(0)));

namespace {
int trailingZeroCount(unsigned int x) {
#if defined(_MSC_VER)
    unsigned long index = 0;
    _BitScanForward(&index, x);
    return static_cast<int>(index);
#else
    return __builtin_ctz(x);
#endif
}

int popcountByte(uint8_t x) {
#if defined(_MSC_VER)
    return __popcnt(static_cast<unsigned int>(x));
#else
    return __builtin_popcount(x);
#endif
}

int evaluateOutputFromKey(const std::vector<uint8_t> &key_bytes, const std::vector<int> &permuted_indices, int a, int b) {
    std::vector<int> register_bits;
    register_bits.reserve(permuted_indices.size());

    for (int idx : permuted_indices) {
        const int byte_idx = idx / 8;
        const int bit_idx = idx % 8;
        const uint8_t byte = key_bytes[byte_idx];
        const int bit = (byte >> (7 - bit_idx)) & 1;
        register_bits.push_back(bit);
    }

    std::vector<int> xor_bits(register_bits.begin(), register_bits.begin() + a);
    std::vector<int> maj_bits(register_bits.begin() + a, register_bits.begin() + a + b);
    return xorMaj(xor_bits, maj_bits);
}
} // namespace

std::string bytesToBitstring(const std::vector<uint8_t> &byte_data) {
    std::string bitstring;
    for (auto byte : byte_data) {
        std::bitset<8> bits(byte);
        bitstring += bits.to_string();
    }
    return bitstring;
}

std::vector<uint8_t> bitstringToBytes(const std::string &bitstring) {
    if (bitstring.empty() || (bitstring.size() % 8 != 0)) {
        throw std::invalid_argument("bitstring length must be a non-zero multiple of 8");
    }

    std::vector<uint8_t> bytes(bitstring.size() / 8, 0);
    for (size_t i = 0; i < bitstring.size(); ++i) {
        const char c = bitstring[i];
        if (c != '0' && c != '1') {
            throw std::invalid_argument("bitstring must contain only '0' or '1'");
        }
        if (c == '1') {
            bytes[i / 8] |= static_cast<uint8_t>(1u << (7 - (i % 8)));
        }
    }
    return bytes;
}

std::vector<std::vector<int>> generateDeterministicSubsets(int num_samples, int N, int n, uint32_t seed) {
    if (num_samples <= 0 || n <= 0 || N <= 0 || n > N) {
        throw std::invalid_argument("invalid subset generation parameters");
    }

    std::vector<std::vector<int>> subsets;
    subsets.reserve(static_cast<size_t>(num_samples));

    std::vector<int> all_indices(N);
    std::iota(all_indices.begin(), all_indices.end(), 0);
    std::mt19937 local_rng(seed);

    for (int i = 0; i < num_samples; ++i) {
        std::vector<int> permuted = all_indices;
        std::shuffle(permuted.begin(), permuted.end(), local_rng);
        permuted.resize(n);
        subsets.push_back(std::move(permuted));
    }

    return subsets;
}

std::vector<uint8_t> generateOutputsFromKey(const std::vector<uint8_t> &key_bytes, const std::vector<std::vector<int>> &subsets, int a, int b) {
    std::vector<uint8_t> outputs;
    outputs.reserve(subsets.size());

    for (const auto &subset : subsets) {
        outputs.push_back(static_cast<uint8_t>(evaluateOutputFromKey(key_bytes, subset, a, b)));
    }

    return outputs;
}

bool verifyCandidateByPRGOutputs(const std::string &candidate_key_bits, const std::vector<std::vector<int>> &subsets, const std::vector<uint8_t> &expected_outputs, int a, int b) {
    if (subsets.size() != expected_outputs.size()) {
        throw std::invalid_argument("subsets and expected_outputs size mismatch");
    }

    const std::vector<uint8_t> candidate_key_bytes = bitstringToBytes(candidate_key_bits);
    const std::vector<uint8_t> candidate_outputs = generateOutputsFromKey(candidate_key_bytes, subsets, a, b);
    return candidate_outputs == expected_outputs;
}

int xorMaj(const std::vector<int> &xor_bits, const std::vector<int> &maj_bits) {
    int xor_result = 0;
    for (int bit : xor_bits) {
        xor_result ^= bit;
    }

    int maj_result = (std::count(maj_bits.begin(), maj_bits.end(), 1) > 31) ? 1 : 0;
    return xor_result ^ maj_result;
}

std::vector<int> generatePermutedIndices(const std::vector<int> &t_indices, int total_size, int total_range) {
    std::vector<int> all_indices(total_range);
    std::iota(all_indices.begin(), all_indices.end(), 0);

    std::vector<int> unique_t_indices = t_indices;
    std::sort(unique_t_indices.begin(), unique_t_indices.end());
    unique_t_indices.erase(std::unique(unique_t_indices.begin(), unique_t_indices.end()), unique_t_indices.end());

    if (unique_t_indices.size() > 64) {
        throw std::invalid_argument("t_indices should have at most 64 unique elements");
    }

    std::vector<int> remaining_indices;
    for (int i = 0; i < total_range; ++i) {
        if (std::find(unique_t_indices.begin(), unique_t_indices.end(), i) == unique_t_indices.end()) {
            remaining_indices.push_back(i);
        }
    }

    std::shuffle(remaining_indices.begin(), remaining_indices.end(), generator);

    std::vector<int> permuted_indices(remaining_indices.begin(), remaining_indices.begin() + total_size - unique_t_indices.size());
    permuted_indices.insert(permuted_indices.end(), unique_t_indices.begin(), unique_t_indices.end());

    std::shuffle(permuted_indices.begin() + 10, permuted_indices.begin() + 74, generator);

    return permuted_indices;
}

std::string recover_key(const std::vector<int> &solution, const std::vector<int> &t_indices, int total_key_size) {
    std::vector<int> key(total_key_size, 0);
    std::vector<int> all_indices(total_key_size);
    std::iota(all_indices.begin(), all_indices.end(), 0);

    std::vector<int> adjusted_t_indices = t_indices;
    std::vector<int> remaining_indices;

    for (int idx : all_indices) {
        if (std::find(adjusted_t_indices.begin(), adjusted_t_indices.end(), idx) == adjusted_t_indices.end()) {
            remaining_indices.push_back(idx);
        }
    }

    for (size_t j = 0; j < remaining_indices.size(); ++j) {
        int i = remaining_indices[j];
        key[i] = solution[j];
    }

    for (int t : adjusted_t_indices) {
        key[t] = 1;
    }

    std::string bitstring;
    for (int bit : key) {
        bitstring += std::to_string(bit);
    }

    return bitstring;
}

void printMatrix(const std::vector<std::vector<int>> &matrix) {
    for (const auto &row : matrix) {
        for (const auto &element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

std::pair<std::vector<std::bitset<MAX_SIZE>>, std::vector<int>> generateMatrixFromInstances(
    const std::vector<std::pair<int, std::vector<int>>> &successfulInstances, const std::vector<int> &tIndices, int totalKeySize) {
    
    std::vector<std::pair<int, std::vector<int>>> selectedInstances = successfulInstances;
    size_t maxRows = 245;
    static std::mt19937 matrix_sampler(static_cast<unsigned int>(std::time(nullptr)));

    if (successfulInstances.size() > maxRows) {
        std::shuffle(selectedInstances.begin(), selectedInstances.end(), matrix_sampler);
        selectedInstances.resize(maxRows);
    }

    std::vector<int> allIndices(totalKeySize);
    std::iota(allIndices.begin(), allIndices.end(), 0);
    std::vector<int> remainingIndices;

    for (int idx : allIndices) {
        if (std::find(tIndices.begin(), tIndices.end(), idx) == tIndices.end()) {
            remainingIndices.push_back(idx);
        }
    }

    size_t rows = selectedInstances.size();
    size_t cols = remainingIndices.size();
    std::vector<std::bitset<MAX_SIZE>> matrix(rows);
    std::vector<int> resultBitsFlipped(rows, 0);

    for (size_t row = 0; row < rows; ++row) {
        auto [resultBit, permutedIndices] = selectedInstances[row];
        std::vector<int> xorIndices(permutedIndices.begin(), permutedIndices.begin() + 10);

        for (int xorIdx : xorIndices) {
            auto it = std::find(remainingIndices.begin(), remainingIndices.end(), xorIdx);
            if (it != remainingIndices.end()) {
                int col = std::distance(remainingIndices.begin(), it);
                if (col >= 0) {
                    matrix[row].set(col);
                }
            }
        }

        resultBitsFlipped[row] = 1 - resultBit;
    }

    return std::make_pair(matrix, resultBitsFlipped);
}

std::vector<std::vector<int>> verifySolutions(const std::vector<std::bitset<MAX_SIZE>> &A, const std::vector<int> &b, const std::vector<std::vector<int>> &solutions) {
    std::vector<std::vector<int>> correct_solutions;

    for (const auto &solution : solutions) {
        bool is_correct = true;
        std::bitset<MAX_SIZE> solution_bitset;

        for (int j = 0; j < MAX_SIZE; ++j) {
            solution_bitset[j] = (j < static_cast<int>(solution.size())) ? solution[j] : 0;
        }

        for (int i = 0; i < A.size(); ++i) {
            int sum = (A[i] & solution_bitset).count() % 2;

            if (sum != b[i]) {
                is_correct = false;
                break;
            }
        }

        if (is_correct) {
            correct_solutions.emplace_back(solution);
        }
    }

    return correct_solutions;
}

bool checkMajPart(const std::vector<int> &guess_indices, const std::vector<int> &permuted_indices) {
    for (int idx : guess_indices) {
        if (std::find(permuted_indices.end() - 64, permuted_indices.end(), idx) == permuted_indices.end()) {
            return false;
        }
    }
    return true;
}

std::vector<std::vector<int>> gaussianEliminationMod2WithFreeVariables(std::vector<std::bitset<MAX_SIZE>> &A, std::vector<int> &b, int rows, int cols) {
    std::vector<int> x(cols, 0);
    std::vector<int> free_variables;
    std::vector<int> leading_variables;

    for (int i = 0; i < std::min(rows, cols); ++i) {
        int pivot_row = -1;
        for (int j = i; j < rows; ++j) {
            if (A[j][i]) {
                pivot_row = j;
                break;
            }
        }

        if (pivot_row == -1) {
            free_variables.push_back(i);
            continue;
        }

        if (pivot_row != i) {
            std::swap(A[i], A[pivot_row]);
            std::swap(b[i], b[pivot_row]);
        }

        leading_variables.push_back(i);

#pragma omp parallel for
        for (int j = i + 1; j < rows; ++j) {
            if (A[j][i]) {
                A[j] ^= A[i];
                b[j] ^= b[i];
            }
        }
    }

    for (int var : free_variables) {
        x[var] = 0;
    }

    for (int i = leading_variables.size() - 1; i >= 0; --i) {
        int var = leading_variables[i];
        int sum = b[var];
        for (int j = var + 1; j < cols; ++j) {
            sum = (sum + A[var][j] * x[j]) % 2;
        }
        x[var] = sum;
    }

    if (free_variables.empty()) {
        return {x};
    }

    std::vector<std::vector<int>> all_solutions;
    int num_free_vars = free_variables.size();
    std::vector<std::vector<int>> dp(1 << num_free_vars, std::vector<int>(cols, 0));
    dp[0] = x;
    all_solutions.push_back(dp[0]);

    for (int combination = 1; combination < (1 << num_free_vars); ++combination) {
        int idx = trailingZeroCount(static_cast<unsigned int>(combination));
        dp[combination] = dp[combination ^ (1 << idx)];

        int var = free_variables[idx];
        dp[combination][var] ^= 1;

        for (int i = leading_variables.size() - 1; i >= 0; --i) {
            int lead_var = leading_variables[i];
            int sum = b[lead_var];
            for (int j = lead_var + 1; j < cols; ++j) {
                sum = (sum + A[lead_var][j] * dp[combination][j]) % 2;
            }
            dp[combination][lead_var] = sum;
        }

        all_solutions.push_back(dp[combination]);
    }

    return all_solutions;
}

std::vector<int> selectRandomIndices(int t) {
    if (t > 256) {
        throw std::invalid_argument("t must be less than or equal to 256");
    }

    std::vector<int> all_indices(256);
    std::iota(all_indices.begin(), all_indices.end(), 0);
    std::shuffle(all_indices.begin(), all_indices.end(), generator);
    return std::vector<int>(all_indices.begin(), all_indices.begin() + t);
}

int hammingWeight(const uint8_t *a, const uint8_t *b, size_t size) {
    int weight = 0;
    for (size_t i = 0; i < size; ++i) {
        weight += popcountByte(static_cast<uint8_t>(a[i] ^ b[i]));
    }
    return weight;
}

void writeResultToFile(std::ofstream &outFile, uint8_t result_bit, const uint8_t *xor_vector, const uint8_t *maj_vector, size_t size) {
    outFile << static_cast<int>(result_bit) << " ";

    for (size_t i = 0; i < size * 8; ++i) {
        if (xor_vector[i / 8] & (1 << (i % 8))) {
            outFile << i << " ";
        }
    }

    for (size_t i = 0; i < size * 8; ++i) {
        if (maj_vector[i / 8] & (1 << (i % 8))) {
            outFile << i << " ";
        }
    }

    outFile << "\n";
}

void processSingleFileWithMemoryMapping(const std::string &inputFilePath, const std::vector<int> &t_indices, std::ofstream &outFile) {
    std::wstring wideInputFilePath = std::wstring(inputFilePath.begin(), inputFilePath.end());

    HANDLE hFile = CreateFileW(wideInputFilePath.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        std::cerr << "cannot open input file" << std::endl;
        return;
    }

    HANDLE hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (hMapFile == NULL) {
        std::cerr << "file mapping failed" << std::endl;
        CloseHandle(hFile);
        return;
    }

    uint8_t *fileData = (uint8_t *)MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);
    if (fileData == NULL) {
        std::cerr << "map view failed" << std::endl;
        CloseHandle(hMapFile);
        CloseHandle(hFile);
        return;
    }

    uint8_t t_indices_vector[32] = {0};
    for (int idx : t_indices) {
        t_indices_vector[idx / 8] |= (1 << (idx % 8));
    }

    size_t offset = 0;
    size_t fileSize = GetFileSize(hFile, NULL);
    while (offset < fileSize) {
        uint8_t result_bit = fileData[offset];
        offset += 1;

        const uint8_t *xor_vector = fileData + offset;
        offset += 32;
        const uint8_t *maj_vector = fileData + offset;
        offset += 32;

        int weight = hammingWeight(maj_vector, t_indices_vector);
        if (weight == 52) {
            writeResultToFile(outFile, result_bit, xor_vector, maj_vector);
        }
    }

    UnmapViewOfFile(fileData);
    CloseHandle(hMapFile);
    CloseHandle(hFile);
}

void processFilesWithMemoryMapping(const std::string &folderPath, const std::vector<int> &t_indices, const std::string &outputFilePath, int num_files) {
    std::ofstream outFile(outputFilePath);
    if (!outFile.is_open()) {
        std::cerr << "cannot open output file" << std::endl;
        return;
    }

    outFile << "t_indices: ";
    for (int idx : t_indices) {
        outFile << idx << " ";
    }
    outFile << "\n";

    for (int i = 1; i <= num_files; ++i) {
        std::string inputFilePath = folderPath + "/output_" + std::to_string(i) + ".bin";
        std::cout << "processing file: " << inputFilePath << std::endl;
        processSingleFileWithMemoryMapping(inputFilePath, t_indices, outFile);
    }

    outFile.close();
}

