#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <iterator>
#include <string>
#include <sstream>  // 用于 istringstream
#include <algorithm>  // 用于 std::count
#include <ctime>

#include "PRNG.h"
#include "verify.h"

using namespace std;

int main()
{
    int a = 10;
    int b = 64;
    vector<uint8_t> key = {
        0x28, 0x33, 0xbe, 0xad, 0x04, 0x9d, 0xd7, 0xb6, 0x34, 0xcf, 0x88, 0x9c, 0x16, 0x5e,
        0x95, 0x57, 0xeb, 0xe3, 0xe8, 0x2d, 0x2a, 0xeb, 0x1f, 0x72, 0xf2, 0x29, 0xc7, 0x8b,
        0x26, 0x5b, 0x99, 0xaf};

    PRNG prng(key, 256, 74);
    string key_bitstring = bytesToBitstring(key);
    int num_1 = 0;
    clock_t start_time = clock();

    while (true)
    {
        vector<int> t_indices = selectRandomIndices(12);

        int success_count = 0;
        vector<pair<int, vector<int>>> correct_instances;
        vector<pair<int, vector<int>>> all_instances;
        int num_instances = 257;

        // 处理文件
        processFilesWithMemoryMapping("data", t_indices, "filtered_output.txt", 500); // 假设读取 500 个文件

        ifstream file("filtered_output.txt");
        if (!file.is_open())
        {
            cerr << "can not open filtered_output.txt" << endl;
            return 1;
        }

        string line;
        getline(file, line);
        cout << "t_indices: " << line << endl;

        while (getline(file, line))
        {
            istringstream iss(line);
            vector<int> parts((std::istream_iterator<int>(iss)), std::istream_iterator<int>());


            if (parts.empty())
                continue;

            int result_bit = parts[0];
            vector<int> permuted_indices(parts.begin() + 1, parts.end());

            // 计算寄存器中的比特
            vector<int> register_bits;
            for (int idx : permuted_indices)
            {
                int byte_idx = idx / 8;
                int bit_idx = idx % 8;
                uint8_t byte = prng.initRegister()[byte_idx];
                int bit = (byte >> (7 - bit_idx)) & 1;
                register_bits.push_back(bit);
            }

            // XOR 部分
            vector<int> xor_bits(register_bits.begin(), register_bits.begin() + a);
            vector<int> maj_bits(register_bits.begin() + a, register_bits.begin() + a + b);

            int xor_result = 0;
            for (int bit : xor_bits)
            {
                xor_result ^= bit;
            }

            int maj_result = (count(maj_bits.begin(), maj_bits.end(), 1) >= maj_bits.size() / 2) ? 1 : 0;
            all_instances.emplace_back(result_bit, permuted_indices);

            if (maj_result == 1)
            {
                success_count++;
                correct_instances.emplace_back(result_bit, permuted_indices);
            }
        }
        file.close();

        int num = 0;
        cout << "Number of instances in correct_instances: " << correct_instances.size() << endl;

        while (num < 309904)
        {
            // 生成矩阵和结果向量
            auto [matrix, results_bits] = generateMatrixFromInstances(all_instances, t_indices, 256);

            vector<vector<int>> solutions = gaussianEliminationMod2WithFreeVariables(matrix, results_bits, 244, 244);

            // 验证解的正确性，将所有有效解保存下来
            vector<vector<int>> valid_solutions = verifySolutions(matrix, results_bits, solutions);
            for (int i = 0; i < valid_solutions.size(); ++i)
            {
                for (int j = 0; j < solutions[i].size(); ++j)
                {
                    valid_solutions[i][j] = 1 - solutions[i][j]; // 将1变成0，0变成1
                }
            }

            // 遍历有效解并进行密钥比较
            for (const auto &valid_sol : valid_solutions)
            {
                string guess_key = recover_key(valid_sol, t_indices);
                if (guess_key == key_bitstring)
                {
                    cout << "The recovered key matches the original key." << endl;
                    cout << "tired " << num << " times" << endl;

                    clock_t end_time = clock();
                    double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;
                    cout << "Elapsed time: " << elapsed_time << " seconds" << endl;
                    return 0;
                }
            }

            num++;
        }

        num_1++;
        cout << "tried " << num_1 << " groups" << endl;
    }

    clock_t end_time = clock();
    double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Elapsed time: " << elapsed_time << " seconds" << endl;
    cout << "failed" << endl;

    return 0;
}
