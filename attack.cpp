#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "prng.h"
#include "verify.h"

int main() {
    const int a = 10;
    const int b = 64;
    const int max_group_trials = 309904;
    const int verify_bits = 4096;

    const std::vector<uint8_t> key = {
        0xb7, 0x24, 0x12, 0xe7, 0x91, 0x57, 0xcf, 0x37, 0xcb, 0xfb, 0xe7, 0xee, 0x7b, 0x5b,
        0x02, 0xd9, 0xcf, 0x50, 0x96, 0x25, 0x3e, 0xfa, 0x99, 0xe4, 0x86, 0x3e, 0x78, 0xd8,
        0x05, 0x5c, 0xd6, 0xf5
    };

    PRNG prng(key, 256, 74);
    const std::vector<uint8_t> register_data = prng.initRegister();

    // Precompute reference outputs once for output-based key verification.
    const auto verify_subsets = generateDeterministicSubsets(verify_bits, 256, 74, 20260228u);
    const auto expected_outputs = generateOutputsFromKey(key, verify_subsets, a, b);

    int tried_groups = 0;
    const clock_t start_time = clock();

    while (true) {
        std::vector<int> t_indices = selectRandomIndices(12);
        std::vector<std::pair<int, std::vector<int>>> all_instances;
        std::vector<std::pair<int, std::vector<int>>> correct_instances;

        processFilesWithMemoryMapping("data", t_indices, "filtered_output.txt", 500);

        std::ifstream file("filtered_output.txt");
        if (!file.is_open()) {
            std::cerr << "cannot open filtered_output.txt" << std::endl;
            return 1;
        }

        std::string line;
        std::getline(file, line);
        std::cout << "t_indices: " << line << std::endl;

        // Keep only equations whose MAJ part is 1.
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<int> parts((std::istream_iterator<int>(iss)), std::istream_iterator<int>());
            if (parts.empty()) {
                continue;
            }

            const int result_bit = parts[0];
            std::vector<int> permuted_indices(parts.begin() + 1, parts.end());

            std::vector<int> register_bits;
            register_bits.reserve(permuted_indices.size());
            for (int idx : permuted_indices) {
                const int byte_idx = idx / 8;
                const int bit_idx = idx % 8;
                const uint8_t byte = register_data[byte_idx];
                const int bit = (byte >> (7 - bit_idx)) & 1;
                register_bits.push_back(bit);
            }

            std::vector<int> maj_bits(register_bits.begin() + a, register_bits.begin() + a + b);
            const int maj_result = (std::count(maj_bits.begin(), maj_bits.end(), 1) >= static_cast<int>(maj_bits.size() / 2)) ? 1 : 0;
            all_instances.emplace_back(result_bit, permuted_indices);
            if (maj_result == 1) {
                correct_instances.emplace_back(result_bit, std::move(permuted_indices));
            }
        }
        file.close();

        std::cout << "Number of instances in correct_instances: " << correct_instances.size() << std::endl;
        std::cout << "Number of instances in all_instances: " << all_instances.size() << std::endl;
        if (all_instances.size() < 245) {
            std::cout << "Not enough equations in this group, retrying next group." << std::endl;
            ++tried_groups;
            continue;
        }

        for (int trial = 0; trial < max_group_trials; ++trial) {
            // Keep attack behavior aligned with the original repository: solve on all filtered equations.
            auto [matrix, rhs] = generateMatrixFromInstances(all_instances, t_indices, 256);
            auto solutions = gaussianEliminationMod2WithFreeVariables(matrix, rhs, 244, 244);
            auto valid_solutions = verifySolutions(matrix, rhs, solutions);

            for (auto &sol : valid_solutions) {
                for (size_t j = 0; j < sol.size(); ++j) {
                    sol[j] = 1 - sol[j];
                }

                const std::string guessed_key = recover_key(sol, t_indices);
                if (verifyCandidateByPRGOutputs(guessed_key, verify_subsets, expected_outputs, a, b)) {
                    std::cout << "Recovered key passed PRG-output verification on " << verify_bits << " bits." << std::endl;
                    std::cout << "solved at trial " << trial << std::endl;
                    const double elapsed = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
                    std::cout << "Elapsed time: " << elapsed << " seconds" << std::endl;
                    return 0;
                }
            }
        }

        ++tried_groups;
        std::cout << "tried " << tried_groups << " groups" << std::endl;
    }

    return 0;
}
