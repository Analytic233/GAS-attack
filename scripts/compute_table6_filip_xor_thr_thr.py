from math import comb, floor, log2


def binomial_tail_prob(total_bits, threshold):
    if threshold <= 0:
        return 1.0
    if threshold > total_bits:
        return 0.0
    return sum(comb(total_bits, i) for i in range(threshold, total_bits + 1)) / (2 ** total_bits)


def q_at_least_r(u, r):
    return sum(comb(u, v) for v in range(r, u + 1)) / (2 ** u)


def overlap_prob_two_thr(n, b1, b2, l, u1, u2):
    if u1 < 0 or u2 < 0:
        return 0.0
    if u1 > l or u1 > b1:
        return 0.0
    if u2 > l - u1 or u2 > b2:
        return 0.0
    if b1 - u1 > n - l:
        return 0.0
    if b2 - u2 > n - l - (b1 - u1):
        return 0.0

    numerator = (
        comb(l, u1)
        * comb(l - u1, u2)
        * comb(n - l, b1 - u1)
        * comb(n - l - (b1 - u1), b2 - u2)
    )
    denominator = comb(n, b1) * comb(n - b1, b2)
    return numerator / denominator


TABLE6_ROWS = [
    {"a": 40, "t1": 26, "b1": 52, "t2": 26, "b2": 52, "n": 360, "r1": 11, "r2": 11, "l": 42, "lambda_bits": 80},
    {"a": 100, "t1": 11, "b1": 22, "t2": 11, "b2": 22, "n": 397, "r1": 7, "r2": 7, "l": 40, "lambda_bits": 80},
    {"a": 54, "t1": 19, "b1": 38, "t2": 19, "b2": 38, "n": 6500, "r1": 7, "r2": 7, "l": 109, "lambda_bits": 128},
    {"a": 54, "t1": 22, "b1": 45, "t2": 23, "b2": 45, "n": 3072, "r1": 9, "r2": 10, "l": 111, "lambda_bits": 128},
    {"a": 70, "t1": 30, "b1": 61, "t2": 31, "b2": 61, "n": 841, "r1": 13, "r2": 14, "l": 61, "lambda_bits": 128},
    {"a": 65, "t1": 47, "b1": 95, "t2": 48, "b2": 96, "n": 830, "r1": 16, "r2": 18, "l": 66, "lambda_bits": 128},
    {"a": 70, "t1": 46, "b1": 93, "t2": 47, "b2": 93, "n": 736, "r1": 16, "r2": 18, "l": 62, "lambda_bits": 128},
    {"a": 60, "t1": 26, "b1": 53, "t2": 27, "b2": 53, "n": 1229, "r1": 11, "r2": 13, "l": 74, "lambda_bits": 128},
    {"a": 160, "t1": 24, "b1": 48, "t2": 24, "b2": 48, "n": 913, "r1": 12, "r2": 12, "l": 62, "lambda_bits": 128},
]


def compute_row(row):
    p1 = binomial_tail_prob(row["b1"] - row["r1"], row["t1"] - row["r1"])
    p2 = binomial_tail_prob(row["b2"] - row["r2"], row["t2"] - row["r2"])
    p = (1.0 + (2.0 * p1 - 1.0) * (2.0 * p2 - 1.0)) / 2.0

    coeff_b = 0.0
    coeff_e = 0.0
    max_u1 = min(row["l"], row["b1"])
    for u1 in range(row["r1"], max_u1 + 1):
        max_u2 = min(row["l"] - u1, row["b2"])
        for u2 in range(row["r2"], max_u2 + 1):
            pu = overlap_prob_two_thr(row["n"], row["b1"], row["b2"], row["l"], u1, u2)
            coeff_b += pu
            coeff_e += pu * q_at_least_r(u1, row["r1"]) * q_at_least_r(u2, row["r2"])

    ep_target = row["n"] - row["l"]
    e_target = ep_target / p
    m = e_target / coeff_e
    d = log2(m)
    b_filtered = m * coeff_b
    e = m * coeff_e
    ep = e * p
    teg = row["l"] + log2(b_filtered)
    ts = row["l"] + ep * (-log2(p))
    tbit = ts + 3 * log2(row["n"] - row["l"])
    tenc = tbit - log2(2 * (row["a"] + row["b1"] + row["b2"]))

    return {
        "a": row["a"],
        "t1": row["t1"],
        "b1": row["b1"],
        "t2": row["t2"],
        "b2": row["b2"],
        "n": row["n"],
        "r1": row["r1"],
        "r2": row["r2"],
        "l": row["l"],
        "p": p,
        "D": d,
        "e": e,
        "ep": ep,
        "Teg": teg,
        "Ts": ts,
        "Tbit": tbit,
        "Tenc": tenc,
        "lambda_bits": row["lambda_bits"],
    }


def main():
    print("Table 6 (paper-style precision)")
    print("a\tt1\tb1\tt2\tb2\tn\tr1\tr2\tl\tp\tD\te\tep\tTeg\tTs\tTbit\tTenc\t2^lambda")
    for row in TABLE6_ROWS:
        result = compute_row(row)
        print(
            f'{result["a"]}\t{result["t1"]}\t{result["b1"]}\t{result["t2"]}\t{result["b2"]}\t{result["n"]}\t'
            f'{result["r1"]}\t{result["r2"]}\t{result["l"]}\t{result["p"]:.4f}\t{result["D"]:.1f}\t'
            f'{floor(result["e"])}\t{floor(result["ep"])}\t{result["Teg"]:.1f}\t{result["Ts"]:.1f}\t'
            f'{result["Tbit"]:.1f}\t{result["Tenc"]:.1f}\t{result["lambda_bits"]}'
        )


if __name__ == "__main__":
    main()
