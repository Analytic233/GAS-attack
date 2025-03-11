#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>  // 用于 256 位 01 向量
#include <iomanip>
#include <random>
#include <string>
#include <sstream>  // 用于生成文件名
#include <sys/stat.h>   // 用于创建目录
#include <sys/types.h>  // 用于目录相关操作

#ifdef _WIN32
#include <direct.h>  // Windows 下用于 _mkdir
#else
#include <unistd.h>
#endif

#include "PRNG.h"  // 引入头文件

using namespace std;

// 新的写入函数：将 bitset<256> 按字节写入文件
void writeBitsetToFile(ofstream &outFile, const bitset<256> &bits) {
    for (size_t i = 0; i < 256; i += 8) {
        uint8_t byte = 0;
        for (int bit = 0; bit < 8; ++bit) {
            if (bits[i + bit]) {
                byte |= (1 << bit); // 设置对应的位
            }
        }
        outFile.write(reinterpret_cast<const char *>(&byte), sizeof(byte));
    }
}

// 检查目录是否存在，如果不存在则创建
void createDataDirectory(const string &directory) {
    struct stat info;
    if (stat(directory.c_str(), &info) != 0) {
        // 如果目录不存在，则创建目录
#ifdef _WIN32
        _mkdir(directory.c_str());  // Windows下创建目录
#else
        mkdir(directory.c_str(), 0777);  // Linux下创建目录
#endif
        std::cout << "create: " << directory << std::endl;
    }
}

int main() {
    // 初始化密钥
    vector<uint8_t> key = {
        0x28, 0x33, 0xbe, 0xad, 0x04, 0x9d, 0xd7, 0xb6, 0x34, 0xcf, 0x88, 0x9c, 0x16, 0x5e,
        0x95, 0x57, 0xeb, 0xe3, 0xe8, 0x2d, 0x2a, 0xeb, 0x1f, 0x72, 0xf2, 0x29, 0xc7, 0x8b,
        0x26, 0x5b, 0x99, 0xaf
    };

    // 设置参数
    int N = 256;                       // N 范围
    int n = 74;                        // 子集大小
    int a = 10;                        // XOR 操作的比特长度
    int b = 64;                        // MAJ 操作的比特长度
    long long num_instances = 10000000000; // 总实例数量

    // 初始化 PRNG 对象
    PRNG prng(key, N, n);

    // 创建 data 文件夹
    string directory = "data";
    createDataDirectory(directory);

    clock_t start_time = clock();

    // 每20000000个实例保存到一个新文件
    int instances_per_file = 20000000;
    int file_index = 1;

    // 打开第一个文件
    stringstream filename;
    filename << directory << "/output_" << file_index << ".bin";  // 添加 .bin 后缀
    ofstream outFile(filename.str(), ios::binary);
    if (!outFile.is_open()) {
        cerr << "open failed!" << endl;
        return 1;
    }

    // 生成实例并写入文件
    for (long long i = 0; i < num_instances; ++i) {
        if (i > 0 && i % instances_per_file == 0) {
            // 关闭当前文件，打开下一个文件
            outFile.close();
            file_index++;
            filename.str("");  // 清空文件名
            filename << directory << "/output_" << file_index << ".bin";  // 添加 .bin 后缀
            outFile.open(filename.str(), ios::binary);
            if (!outFile.is_open()) {
                cerr << "open failed!" << endl;
                return 1;
            }
        }

        // 生成伪随机密钥流
        auto [bit, permuted_indices] = prng.generateKeystream(a, b);

        // 初始化两个 256 位的 01 向量
        bitset<256> xor_vector;
        bitset<256> maj_vector;

        // 将前 a 个索引设置为 XOR 的 1 位
        for (int j = 0; j < a; ++j) {
            xor_vector.set(permuted_indices[j], 1); // 将指定索引设置为1
        }

        // 将接下来 b 个索引设置为 MAJ 的 1 位
        for (int j = a; j < a + b; ++j) {
            maj_vector.set(permuted_indices[j], 1); // 将指定索引设置为1
        }

        // 写入 result_bit (用一个字节表示 0 或 1)
        uint8_t result_bit = static_cast<uint8_t>(bit);
        outFile.write(reinterpret_cast<const char *>(&result_bit), sizeof(result_bit));

        // 写入 XOR 256 位 01 向量（以二进制形式）
        writeBitsetToFile(outFile, xor_vector);

        // 写入 MAJ 256 位 01 向量（以二进制形式）
        writeBitsetToFile(outFile, maj_vector);

        // 可选：添加进度显示，每生成 100 万个实例显示一次进度
        if (i % 1000000 == 0) {
            cout << "created " << i << " outputs" << endl;
        }
    }

    // 关闭文件
    outFile.close();

    cout << "created " << num_instances << " outputs" << endl;
    clock_t end_time = clock();
    double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Elapsed time: " << elapsed_time << " seconds" << endl;

    return 0;
}
