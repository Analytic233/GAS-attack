#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <map>

using namespace std;

// --- 1. parameter setting (from Table 4 ) ---
const int N =1777;
const int LEN_A = 70;
const int LEN_B = 122;
const int THR_T = 61;
const int LOC_N = 192;
const int GUESS_L = 61;
const int FIX_R = 21;
const int SAMPLE_SIZE = 10000000;
mt19937 rng(1);

double log_nCr(int n, int k) {
    if (k < 0 || k > n) return -1.0/0.0; // log(0) -> -inf
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}

double hypergeometric_prob(int u) {
    return exp(log_nCr(GUESS_L, u) + log_nCr(N - GUESS_L, LEN_B - u) - log_nCr(N, LEN_B));
}

class FiLIP_Oracle {
    vector<uint8_t> key;
public:
    FiLIP_Oracle() { key.resize(N); for(auto& bit : key) bit = rng()%2; }
    uint8_t get_key_bit(int idx) { return key[idx]; }
    struct Sample { int y; vector<int> indices; vector<uint8_t> noise; };

    Sample generate_sample_with_intersection(const vector<int>& l_indices, int u) {
        Sample s; s.indices.resize(LOC_N); s.noise.resize(LOC_N);
        vector<int> l_shuffled = l_indices; shuffle(l_shuffled.begin(), l_shuffled.end(), rng);
        
        vector<int> non_l_pool; vector<bool> in_L(N, false); for(int idx : l_indices) in_L[idx] = true;
        for(int i=0; i<N; ++i) if(!in_L[i]) non_l_pool.push_back(i);
        shuffle(non_l_pool.begin(), non_l_pool.end(), rng);

        vector<int> thr_indices;
        for(int i=0; i<u; ++i) thr_indices.push_back(l_shuffled[i]);
        for(int i=0; i<LEN_B - u; ++i) thr_indices.push_back(non_l_pool[i]);
        shuffle(thr_indices.begin(), thr_indices.end(), rng);
        
        vector<int> xor_indices;
        for(int i=LEN_B - u; i < (LEN_B - u) + LEN_A; ++i) xor_indices.push_back(non_l_pool[i]);

        for(int i=0; i<LEN_A; ++i) s.indices[i] = xor_indices[i];
        for(int i=0; i<LEN_B; ++i) s.indices[LEN_A + i] = thr_indices[i];

        for(int i=0; i<LOC_N; ++i) s.noise[i] = rng() % 2;
        int xor_val=0; for(int i=0; i<LEN_A; ++i) xor_val ^= (key[s.indices[i]] ^ s.noise[i]);
        int thr_sum=0; for(int i=LEN_A; i<LOC_N; ++i) thr_sum += (key[s.indices[i]] ^ s.noise[i]);
        s.y = xor_val ^ (thr_sum >= THR_T ? 1 : 0);
        return s;
    }
};


int main() {
    cout << "=== FiLIP High-Fidelity Statistical Verification (Corrected) ===" << endl;
    
    vector<double> conditional_probs;
    vector<int> u_values;
    double total_prob_in_range = 0;
    for (int u = FIX_R; u <= GUESS_L; ++u) {
        double p = hypergeometric_prob(u);
        conditional_probs.push_back(p);
        total_prob_in_range += p;
        u_values.push_back(u);
    }
    
    for(auto& p : conditional_probs) p /= total_prob_in_range;
    std::discrete_distribution<> dist(conditional_probs.begin(), conditional_probs.end());

    // Initialize Oracle
    FiLIP_Oracle oracle;
    vector<int> L_indices(N); iota(L_indices.begin(), L_indices.end(), 0);
    shuffle(L_indices.begin(), L_indices.end(), rng); L_indices.resize(GUESS_L);
    vector<uint8_t> good_key(N); for(int i=0; i<N; ++i) good_key[i]=oracle.get_key_bit(i);
    vector<uint8_t> bad_key = good_key; for(int idx : L_indices) if(idx%2==0) bad_key[idx]^=1;
    vector<bool> in_L(N, false); for(int idx : L_indices) in_L[idx] = true;
    
    int good_matches=0, good_total=0, bad_matches=0, bad_total=0;
    cout << "\n[Running] Analyzing " << SAMPLE_SIZE << " high-fidelity samples..." << endl;
    
    for(int i=0; i<SAMPLE_SIZE; ++i) {
        int u = u_values[dist(rng)];
        auto s = oracle.generate_sample_with_intersection(L_indices, u);

        // --- Good Group ---
        int valid_good=0;
        for(int k=LEN_A; k<LOC_N; ++k) if(in_L[s.indices[k]] && (good_key[s.indices[k]]^s.noise[k])==1) valid_good++;
        if(valid_good >= FIX_R) {
            good_total++;
            int lhs=0; for(int k=0; k<LEN_A; ++k) lhs^=(good_key[s.indices[k]]^s.noise[k]);
            if((s.y^lhs)==1) good_matches++;
        }
        
        // --- Bad Group ---
        int valid_bad=0;
        for(int k=LEN_A; k<LOC_N; ++k) if(in_L[s.indices[k]] && (bad_key[s.indices[k]]^s.noise[k])==1) valid_bad++;
        if(valid_bad >= FIX_R) {
            bad_total++;
            int lhs=0; for(int k=0; k<LEN_A; ++k) lhs^=(good_key[s.indices[k]]^s.noise[k]);
            if((s.y^lhs)==1) bad_matches++;
        }
    }

    // output
    double good_p = (good_total>0) ? (double)good_matches/good_total : 0.0;
    double bad_p = (bad_total>0) ? (double)bad_matches/bad_total : 0.0;
    
    cout << "\n--- Result ---" << endl;
    cout << "Total Samples Generated (filtered): " << SAMPLE_SIZE << endl;
    cout << "------------------------------------------" << endl;
    cout << setw(20) << "Group Type" << " | " << setw(10) << "Equations" << " | " << "Probability" << endl;
    cout << "------------------------------------------" << endl;
    cout << setw(20) << "Good Group" << " | " << setw(10) << good_total << " | " << fixed << setprecision(5) << good_p << endl;
    cout << setw(20) << "Bad Group" << " | " << setw(10) << bad_total << " | " << bad_p << endl;
    cout << "------------------------------------------" << endl;

    return 0;
}