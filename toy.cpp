#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <bitset>
#include <ctime>
#include <iomanip>

// --- Configuration ---
const int N = 60;          // Key size (reduced from 256)
const int LOC_N = 20;      // Locality n (reduced from 74)
const int LEN_A = 5;       // XOR part length
const int LEN_B = 15;      // MAJ part length
const int FIX_R = 5;       // Number of bits to fix/guess
const int EXPECTED_EQ = 55; // Equations needed to solve (approx 2*N)
const int SAMPLE_SIZE = 150000; // Total samples to generate

using namespace std;

// --- PRNG Implementation ---
class PRNG {
    vector<uint8_t> key;
    mt19937 gen;
    int n_size;

public:
    PRNG(const vector<uint8_t>& k, int n) : key(k), n_size(n), gen(43) {}

    pair<int, vector<int>> generate_instance() {
        vector<int> indices(key.size());
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), gen);
        
        vector<int> subset;
        subset.reserve(n_size);
        for(int i=0; i<n_size; ++i) subset.push_back(indices[i]);

        // Extract bits
        vector<int> bits;
        for(int idx : subset) bits.push_back(key[idx]);

        // XOR part
        int xor_val = 0;
        for(int i=0; i<LEN_A; ++i) xor_val ^= bits[i];

        // MAJ part
        int maj_count = 0;
        for(int i=LEN_A; i<LEN_A+LEN_B; ++i) maj_count += bits[i];
        int maj_val = (maj_count > LEN_B/2) ? 1 : 0;

        return {xor_val ^ maj_val, subset};
    }
};

// --- Gaussian Elimination ---
// Returns valid solutions for x
vector<vector<int>> solve_system(vector<bitset<N>>& A, vector<int>& b, int rows) {
    int cols = N;
    vector<int> x(cols, 0);
    vector<int> pivot_row(cols, -1);
    
    // Forward elimination
    int p_row = 0;
    for(int col=0; col<cols && p_row < rows; ++col) {
        int sel = -1;
        for(int r=p_row; r<rows; ++r) {
            if(A[r][col]) { sel = r; break; }
        }
        if(sel == -1) continue;

        swap(A[p_row], A[sel]);
        swap(b[p_row], b[sel]);
        pivot_row[col] = p_row;

        for(int r=0; r<rows; ++r) {
            if(r != p_row && A[r][col]) {
                A[r] ^= A[p_row];
                b[r] ^= b[p_row];
            }
        }
        p_row++;
    }

    for(int i=cols-1; i>=0; --i) {
        if(pivot_row[i] != -1) {
            int val = b[pivot_row[i]];
            for(int j=i+1; j<cols; ++j) {
                if(A[pivot_row[i]][j]) val ^= x[j];
            }
            x[i] = val;
        }
    }
    
    // Simple verification check inside solver to return strictly valid
    return {x}; 
}

// --- Verification ---
bool verify_key(const vector<int>& candidate, const vector<uint8_t>& original) {
    for(int i=0; i<N; ++i) {
        if(candidate[i] != original[i]) return false;
    }
    return true;
}

int main() {
    cout << "=== Toy Example: Attack on Goldreich's PRG ===" << endl;
    cout << "Parameters: N=" << N << ", n=" << LOC_N << ", a=" << LEN_A << ", b=" << LEN_B << endl;
    cout << "Attack Guess Size: r=" << FIX_R << endl;

    // 1. Setup Oracle (The Target)
    vector<uint8_t> secret_key(N);
    mt19937 rng(time(0));
    cout << "Secret Key: ";
    for(int i=0; i<N; ++i) {
        secret_key[i] = rng() % 2;
        cout << (int)secret_key[i];
    }
    cout << endl;

    PRNG oracle(secret_key, LOC_N);

    // 2. Generate Data Pool (In-Memory)
    cout << "[*] Generating " << SAMPLE_SIZE << " samples..." << endl;
    struct Sample { int y; vector<int> indices; };
    vector<Sample> pool;
    pool.reserve(SAMPLE_SIZE);
    for(int i=0; i<SAMPLE_SIZE; ++i) {
        auto res = oracle.generate_instance();
        pool.push_back({res.first, res.second});
    }

    // 3. Attack Loop
    int max_groups = 1000; // Try up to 1000 random groups
    cout << "[*] Starting Group-and-Solve Attack..." << endl;

    for(int g=0; g<max_groups; ++g) {
        // A. Pick random 'r' indices to guess as 1s
        vector<int> guess_indices(N);
        iota(guess_indices.begin(), guess_indices.end(), 0);
        shuffle(guess_indices.begin(), guess_indices.end(), rng);
        guess_indices.resize(FIX_R);

        // B. Filter samples
        vector<Sample> filtered;
        for(const auto& s : pool) {
            // Get MAJ part indices
            vector<int> maj_indices;
            for(int k=LEN_A; k<LEN_A+LEN_B; ++k) maj_indices.push_back(s.indices[k]);

            // Check if guess_indices is a subset of maj_indices
            bool is_subset = true;
            for(int gi : guess_indices) {
                bool found = false;
                for(int mi : maj_indices) {
                    if(mi == gi) { found = true; break; }
                }
                if(!found) { is_subset = false; break; }
            }

            if(is_subset) {
                filtered.push_back(s);
            }
        }

        if(filtered.size() < EXPECTED_EQ) {
            if (g % 100 == 0) cout << "Group " << g << ": Not enough equations (" << filtered.size() << ")" << endl;
            continue; 
        }

        // C. Construct Linear System
        // If guess is correct (all 1s), MAJ outputs 1 with high probability.
        // y = XOR + MAJ. If MAJ -> 1, then XOR = y ^ 1.
        vector<bitset<N>> Matrix;
        vector<int> Rhs;
        
        // Take a subset of equations to solve
        for(int i=0; i<min((int)filtered.size(), 200); ++i) {
            bitset<N> row;
            // Add XOR indices to row
            for(int k=0; k<LEN_A; ++k) {
                row.flip(filtered[i].indices[k]);
            }
            Matrix.push_back(row);
            Rhs.push_back(filtered[i].y ^ 1); // target = y + 1
        }

        // D. Solve
        vector<vector<int>> candidates = solve_system(Matrix, Rhs, Matrix.size());

        // E. Verify
        for(auto& sol : candidates) {
            // Fix the guessed bits to 1
            for(int gi : guess_indices) sol[gi] = 1;
            
            if(verify_key(sol, secret_key)) {
                cout << "\n[SUCCESS] Key Recovered in Group " << g << "!" << endl;
                cout << "Recovered: ";
                for(int b : sol) cout << b;
                cout << endl;
                return 0;
            }
        }
        if (g % 50 == 0) cout << "." << flush;
    }

    cout << "\n[FAIL] Could not recover key in " << max_groups << " attempts." << endl;
    return 1;
}