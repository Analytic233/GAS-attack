#include "prng.h"

#include <algorithm>
#include <ctime>
#include <numeric>

PRNG::PRNG(std::vector<uint8_t> key, int N, int n)
    : key(std::move(key)), N(N), n(n), generator(static_cast<unsigned int>(std::time(nullptr))) {
    register_data = initRegister();
}

std::vector<uint8_t> PRNG::initRegister() {
    return key;
}

std::vector<int> PRNG::generate_subset() {
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), generator);
    indices.resize(n);
    return indices;
}

int PRNG::xorMaj(const std::vector<int>& xorBits, const std::vector<int>& majBits) {
    int xorResult = 0;
    for (int bit : xorBits) {
        xorResult ^= bit;
    }

    const int ones = static_cast<int>(std::count(majBits.begin(), majBits.end(), 1));
    const int majResult = (ones >= static_cast<int>(majBits.size() / 2)) ? 1 : 0;
    return xorResult ^ majResult;
}

std::pair<int, std::vector<int>> PRNG::generateKeystream(int a, int b) {
    std::vector<int> permutedIndices = generate_subset();

    std::vector<int> registerBits;
    registerBits.reserve(permutedIndices.size());
    for (int idx : permutedIndices) {
        const int byteIdx = idx / 8;
        const int bitIdx = idx % 8;
        const uint8_t byte = register_data[byteIdx];
        const int bit = (byte >> (7 - bitIdx)) & 1;
        registerBits.push_back(bit);
    }

    std::vector<int> xorBits(registerBits.begin(), registerBits.begin() + a);
    std::vector<int> majBits(registerBits.begin() + a, registerBits.begin() + a + b);

    const int filteredOutput = xorMaj(xorBits, majBits);
    return {filteredOutput, permutedIndices};
}
