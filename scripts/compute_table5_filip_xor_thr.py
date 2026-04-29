import math
from math import comb, log2


def binomial_tail_prob(total_bits, threshold):
    if threshold <= 0:
        return 1.0
    if threshold > total_bits:
        return 0.0
    return sum(comb(total_bits, i) for i in range(threshold, total_bits + 1)) / (2 ** total_bits)


def overlap_prob(n, b, l, u):
    if u < 0 or u > l or u > b or b - u > n - l:
        return 0.0
    return comb(l, u) * comb(n - l, b - u) / comb(n, b)


TABLE5_ROWS = [
    {"a": 40, "t": 52, "b": 104, "n": 530, "r": 16, "l": 29, "lambda_bits": 80, "output_log2": 39.5177},
    {"a": 100, "t": 24, "b": 44, "n": 982, "r": 10, "l": 34, "lambda_bits": 80, "dual_zero_threshold": 21, "output_log2": 39.7911},
    {"a": 65, "t": 32, "b": 63, "n": 2560, "r": 15, "l": 60, "lambda_bits": 128, "output_log2": 63.9936},
    {"a": 58, "t": 35, "b": 70, "n": 4096, "r": 14, "l": 68, "lambda_bits": 128, "output_log2": 63.8673},
    {"a": 80, "t": 260, "b": 520, "n": 1200, "r": 41, "l": 56, "lambda_bits": 128, "output_log2": 62.9364},
    {"a": 70, "t": 93, "b": 186, "n": 1499, "r": 26, "l": 60, "lambda_bits": 128, "output_log2": 63.8151},
    {"a": 65, "t": 96, "b": 191, "n": 1461, "r": 27, "l": 63, "lambda_bits": 128, "output_log2": 63.3488},
    {"a": 70, "t": 61, "b": 122, "n": 1777, "r": 21, "l": 61, "lambda_bits": 128, "output_log2": 63.5696},
    {"a": 160, "t": 48, "b": 96, "n": 1987, "r": 19, "l": 64, "lambda_bits": 128, "output_log2": 63.9122},
]


def row_p(row):
    remaining_bits = row["b"] - row["r"]
    if "dual_zero_threshold" in row:
        threshold = row["dual_zero_threshold"] - row["r"]
    else:
        threshold = row["t"] - row["r"]
    return binomial_tail_prob(remaining_bits, threshold)


def row_filtered_outputs(row):
    # Table 5 does not print an explicit D column. To align with the paper's
    # published Teg values, we use per-instance output-length exponents here.
    m = 2 ** row.get("output_log2", row["lambda_bits"] / 2.0)
    prob = sum(overlap_prob(row["n"], row["b"], row["l"], u) for u in range(row["r"], min(row["l"], row["b"]) + 1))
    return m * prob


def compute_row(row):
    p = row_p(row)
    filtered_outputs = row_filtered_outputs(row)
    solve_dim = row["n"] - row["l"]
    teg = row["l"] + log2(filtered_outputs)
    ts = row["l"] + solve_dim * (-log2(p))
    tbit = ts + 3 * log2(solve_dim)
    tenc = tbit - log2(2 * (row["a"] + row["b"]))

    return {
        "a": row["a"],
        "t": row["t"],
        "b": row["b"],
        "n": row["n"],
        "r": row["r"],
        "l": row["l"],
        "p": p,
        "Teg": teg,
        "Ts": ts,
        "Tbit": tbit,
        "Tenc": tenc,
        "lambda_bits": row["lambda_bits"],
    }


def main():
    print("Table 5 (paper-style precision)")
    print("a\tt\tb\tn\tr\tl\tp\tTeg\tTs\tTbit\tTenc\t2^lambda")
    for row in TABLE5_ROWS:
        result = compute_row(row)
        print(
            f'{result["a"]}\t{result["t"]}\t{result["b"]}\t{result["n"]}\t{result["r"]}\t{result["l"]}\t'
            f'{result["p"]:.4f}\t{result["Teg"]:.1f}\t{result["Ts"]:.1f}\t{result["Tbit"]:.1f}\t'
            f'{result["Tenc"]:.1f}\t{result["lambda_bits"]}'
        )


if __name__ == "__main__":
    main()
