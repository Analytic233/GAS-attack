#include <bitset>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#include <sys/stat.h>
#include <sys/types.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#include "prng.h"

using namespace std;

// Write bitset<256> as 32 bytes (little-endian by bit index in each byte).
void writeBitsetToFile(ofstream &outFile, const bitset<256> &bits) {
    for (size_t i = 0; i < 256; i += 8) {
        uint8_t byte = 0;
        for (int bit = 0; bit < 8; ++bit) {
            if (bits[i + bit]) {
                byte |= static_cast<uint8_t>(1u << bit);
            }
        }
        outFile.write(reinterpret_cast<const char *>(&byte), sizeof(byte));
    }
}

void createDataDirectory(const string &directory) {
    struct stat info;
    if (stat(directory.c_str(), &info) != 0) {
#ifdef _WIN32
        _mkdir(directory.c_str());
#else
        mkdir(directory.c_str(), 0777);
#endif
        cout << "create: " << directory << endl;
    }
}

int main() {
    vector<uint8_t> key = {
        0xb7, 0x24, 0x12, 0xe7, 0x91, 0x57, 0xcf, 0x37, 0xcb, 0xfb, 0xe7, 0xee, 0x7b, 0x5b,
        0x02, 0xd9, 0xcf, 0x50, 0x96, 0x25, 0x3e, 0xfa, 0x99, 0xe4, 0x86, 0x3e, 0x78, 0xd8,
        0x05, 0x5c, 0xd6, 0xf5
    };

    const int N = 256;
    const int n = 74;
    const int a = 10;
    const int b = 64;
    const long long num_instances = 10000000000LL;
    const int instances_per_file = 20000000;

    PRNG prng(key, N, n);
    const string directory = "data";
    createDataDirectory(directory);

    const clock_t start_time = clock();
    int file_index = 1;

    stringstream filename;
    filename << directory << "/output_" << file_index << ".bin";
    ofstream outFile(filename.str(), ios::binary);
    if (!outFile.is_open()) {
        cerr << "open failed!" << endl;
        return 1;
    }

    for (long long i = 0; i < num_instances; ++i) {
        if (i > 0 && i % instances_per_file == 0) {
            outFile.close();
            ++file_index;
            filename.str("");
            filename.clear();
            filename << directory << "/output_" << file_index << ".bin";
            outFile.open(filename.str(), ios::binary);
            if (!outFile.is_open()) {
                cerr << "open failed!" << endl;
                return 1;
            }
        }

        auto [bit, permuted_indices] = prng.generateKeystream(a, b);

        bitset<256> xor_vector;
        bitset<256> maj_vector;

        for (int j = 0; j < a; ++j) {
            xor_vector.set(permuted_indices[j], 1);
        }
        for (int j = a; j < a + b; ++j) {
            maj_vector.set(permuted_indices[j], 1);
        }

        uint8_t result_bit = static_cast<uint8_t>(bit);
        outFile.write(reinterpret_cast<const char *>(&result_bit), sizeof(result_bit));
        writeBitsetToFile(outFile, xor_vector);
        writeBitsetToFile(outFile, maj_vector);

        if (i % 1000000 == 0) {
            cout << "created " << i << " outputs" << endl;
        }
    }

    outFile.close();
    cout << "created " << num_instances << " outputs" << endl;
    const clock_t end_time = clock();
    const double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Elapsed time: " << elapsed_time << " seconds" << endl;

    return 0;
}
