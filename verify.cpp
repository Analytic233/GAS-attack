#include "verify.h"
#include "PRNG.h"
#include <algorithm>
#include <bitset>
#include <chrono>
#include <random>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <windows.h>

std::mt19937 generator(static_cast<unsigned int>(std::time(0)));

// 将字节转换为比特串
std::string bytesToBitstring(const std::vector<uint8_t> &byte_data) {
    std::string bitstring;
    for (auto byte : byte_data) {
        std::bitset<8> bits(byte);
        bitstring += bits.to_string();
    }
    return bitstring;
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

    std::random_device rd;
    std::mt19937 generator(rd());

    if (successfulInstances.size() > maxRows) {
        std::shuffle(selectedInstances.begin(), selectedInstances.end(), generator);
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
            solution_bitset[j] = solution[j];
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

    for (int combination = 1; combination < (1 << num_free_vars); ++combination) {
        int idx = __builtin_ctz(combination);
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
        weight += __builtin_popcount(a[i] ^ b[i]);
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
        std::cerr << "无法打开输入文件！" << std::endl;
        return;
    }

    HANDLE hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (hMapFile == NULL) {
        std::cerr << "文件映射失败！" << std::endl;
        CloseHandle(hFile);
        return;
    }

    uint8_t *fileData = (uint8_t *)MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);
    if (fileData == NULL) {
        std::cerr << "文件映射视图失败！" << std::endl;
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
        if (weight == 51) {
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
        std::cerr << "无法打开输出文件！" << std::endl;
        return;
    }

    outFile << "t_indices: ";
    for (int idx : t_indices) {
        outFile << idx << " ";
    }
    outFile << "\n";

    for (int i = 1; i <= num_files; ++i) {
        std::string inputFilePath = folderPath + "/output_" + std::to_string(i) + ".bin";
        std::cout << "正在处理文件: " << inputFilePath << std::endl;
        processSingleFileWithMemoryMapping(inputFilePath, t_indices, outFile);
    }

    outFile.close();
}
