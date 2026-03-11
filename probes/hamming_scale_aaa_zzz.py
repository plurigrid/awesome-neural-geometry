"""
world:// hamming scale aaa-zzz probe — trigram Hamming geometry

Explores the Hamming geometry of all 17576 three-letter trigrams (aaa..zzz)
under two encodings:

  Binary:  15 bits (5 per letter, concatenated)
  Ternary:  9 trits (3 per letter, concatenated)

Full 17576x17576 distance matrix is O(300M) — too expensive.
Instead we use analytical formulas, sampling, and row-sum equivalence classes.

Scaling analysis: 1-letter (26) → 2-letter (676) → 3-letter (17576).
"""

import random
import itertools
from collections import Counter

random.seed(42)

# ---------------------------------------------------------------------------
# Encoding helpers
# ---------------------------------------------------------------------------

def to_digits(val: int, base: int, width: int) -> tuple[int, ...]:
    """Convert val to a fixed-width tuple of digits in the given base."""
    digits = []
    for _ in range(width):
        digits.append(val % base)
        val //= base
    return tuple(reversed(digits))


def encode_trigram(tri: str, base: int, width_per_letter: int) -> tuple[int, ...]:
    """Encode a trigram by concatenating per-letter codes."""
    parts = []
    for ch in tri:
        parts.extend(to_digits(ord(ch) - ord('a'), base, width_per_letter))
    return tuple(parts)


def hamming(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    """Hamming distance between two symbol tuples."""
    return sum(x != y for x, y in zip(a, b))


def all_trigrams() -> list[str]:
    """Generate all 26^3 = 17576 trigrams in alphabetical order."""
    return [chr(a) + chr(b) + chr(c)
            for a in range(ord('a'), ord('z') + 1)
            for b in range(ord('a'), ord('z') + 1)
            for c in range(ord('a'), ord('z') + 1)]


# ---------------------------------------------------------------------------
# 1. Analytical mean pairwise Hamming distance
# ---------------------------------------------------------------------------

def analytical_mean(base: int, width_per_letter: int, n_letters: int = 26, n_positions: int = 3):
    """
    For independent uniform letters, the expected Hamming distance between
    two random codewords = sum over all bit/trit positions of P(disagree).

    Each trigram has n_positions letters (3).  Each letter maps to
    width_per_letter digits in the given base.  Letter values are 0..25,
    each equally likely.

    For a single letter position, digit d at bit-position k:
      digit_k(letter) = (letter // base^(width_per_letter-1-k)) % base

    P(two random letters disagree at digit k) =
        1 - sum_v P(digit_k = v)^2

    Total expected Hamming distance = n_positions * sum_k P(disagree at k)
    """
    total_expected = 0.0
    for k in range(width_per_letter):
        # Count how many of the 26 letters have each digit value at position k
        counts = Counter()
        for letter_val in range(n_letters):
            d = to_digits(letter_val, base, width_per_letter)
            counts[d[k]] += 1
        # P(agree) = sum_v (count_v/26)^2
        p_agree = sum((c / n_letters) ** 2 for c in counts.values())
        p_disagree = 1.0 - p_agree
        total_expected += p_disagree

    # Multiply by number of letter positions in the trigram
    return n_positions * total_expected


def per_position_disagree_probs(base: int, width_per_letter: int, n_letters: int = 26):
    """Return P(disagree) for each digit position within one letter."""
    probs = []
    for k in range(width_per_letter):
        counts = Counter()
        for letter_val in range(n_letters):
            d = to_digits(letter_val, base, width_per_letter)
            counts[d[k]] += 1
        p_agree = sum((c / n_letters) ** 2 for c in counts.values())
        probs.append(1.0 - p_agree)
    return probs


# ---------------------------------------------------------------------------
# 2. Sampled verification
# ---------------------------------------------------------------------------

def sampled_mean(trigrams: list[str], base: int, wpl: int, n_pairs: int = 1000):
    """Sample n_pairs random pairs, compute mean Hamming distance."""
    n = len(trigrams)
    total = 0
    dists = []
    for _ in range(n_pairs):
        i, j = random.sample(range(n), 2)
        a = encode_trigram(trigrams[i], base, wpl)
        b = encode_trigram(trigrams[j], base, wpl)
        d = hamming(a, b)
        dists.append(d)
        total += d
    return total / n_pairs, dists


# ---------------------------------------------------------------------------
# 3. Row-sum equivalence classes (analytical via per-letter row sums)
# ---------------------------------------------------------------------------

def letter_row_sums(base: int, wpl: int, n_letters: int = 26):
    """
    For a single letter, compute row-sum = sum of Hamming distances
    from that letter's code to all other 25 letters' codes.
    """
    codes = [to_digits(v, base, wpl) for v in range(n_letters)]
    sums = []
    for i in range(n_letters):
        s = sum(hamming(codes[i], codes[j]) for j in range(n_letters) if j != i)
        sums.append(s)
    return sums


def trigram_row_sum_from_letter_sums(tri: str, letter_sums: list[int],
                                      n_letters: int = 26):
    """
    Row-sum of a trigram = sum of Hamming distances to all 17576 others.

    For trigram (L1, L2, L3), the distance to (M1, M2, M3) is
    d(L1,M1) + d(L2,M2) + d(L3,M3).

    Summing over all M1,M2,M3:
      sum_{M} d(L,M) = N^2 * sum_others d(L1,m) + N^2 * sum_others d(L2,m) + N^2 * sum_others d(L3,m)

    Wait — more carefully:
      sum_{M1,M2,M3} [d(L1,M1) + d(L2,M2) + d(L3,M3)]
        = N^2 * sum_{M1} d(L1,M1) + N^2 * sum_{M2} d(L2,M2) + N^2 * sum_{M3} d(L3,M3)

    where each inner sum is over all 26 letters (including self, giving 0 for self).
    """
    n2 = n_letters ** 2  # 676
    vals = [ord(ch) - ord('a') for ch in tri]
    # letter_sums[v] = sum of distances to all OTHER 25 letters
    # We need sum over all 26 including self (self contributes 0)
    # So letter_sums_full[v] = letter_sums[v] (since d(v,v)=0)
    return n2 * sum(letter_sums[v] for v in vals)


def compute_row_sum_classes(letter_sums: list[int], n_letters: int = 26, n_positions: int = 3):
    """
    Compute equivalence classes of trigrams by row-sum.

    Since row_sum(L1,L2,L3) = N^2 * (ls[L1] + ls[L2] + ls[L3]),
    the classes are determined by the multiset {ls[L1], ls[L2], ls[L3]}.

    Actually, they're determined by the sum ls[L1] + ls[L2] + ls[L3],
    since N^2 is a constant multiplier.
    """
    # Find all distinct sums ls[L1] + ls[L2] + ls[L3]
    ls_values = letter_sums  # indexed by letter value 0..25
    sum_counts = Counter()
    for a in range(n_letters):
        for b in range(n_letters):
            for c in range(n_letters):
                s = ls_values[a] + ls_values[b] + ls_values[c]
                sum_counts[s] += 1

    return sum_counts


# ---------------------------------------------------------------------------
# 4. Fibonacci-indexed subset: path costs
# ---------------------------------------------------------------------------

def fibonacci_indices(limit: int) -> list[int]:
    """Generate Fibonacci numbers (starting 0,1,1,2,...) less than limit."""
    fibs = []
    a, b = 0, 1
    while a < limit:
        fibs.append(a)
        a, b = b, a + b
    return fibs


def path_cost(trigrams: list[str], indices: list[int],
              base: int, wpl: int) -> int:
    """Sum of consecutive Hamming distances along the path."""
    total = 0
    for k in range(len(indices) - 1):
        a = encode_trigram(trigrams[indices[k]], base, wpl)
        b = encode_trigram(trigrams[indices[k + 1]], base, wpl)
        total += hamming(a, b)
    return total


# ---------------------------------------------------------------------------
# 5. Scaling analysis: 1-letter, 2-letter, 3-letter
# ---------------------------------------------------------------------------

def scaling_analysis():
    """How metrics scale from 1-letter to 2-letter to 3-letter."""
    print(f"\n{'='*64}")
    print(f"  SCALING ANALYSIS: 1-letter → 2-letter → 3-letter")
    print(f"{'='*64}\n")

    for base, wpl, label in [(2, 5, "Binary"), (3, 3, "Ternary")]:
        print(f"  --- {label} (base={base}, {wpl} digits/letter) ---\n")
        print(f"  {'N-gram':<8s} {'#items':>7s} {'Width':>6s} "
              f"{'E[d]':>8s} {'E[d]/w':>8s} {'#RowSum classes':>16s}")
        print(f"  {'-'*8} {'-'*7} {'-'*6} {'-'*8} {'-'*8} {'-'*16}")

        ls = letter_row_sums(base, wpl)

        for n_pos in [1, 2, 3]:
            n_items = 26 ** n_pos
            width = wpl * n_pos
            e_d = analytical_mean(base, wpl, n_letters=26, n_positions=n_pos)

            # Row-sum classes: determined by sum of n_pos letter row-sums
            if n_pos == 1:
                class_counts = Counter(Counter(ls).values())
                n_classes = len(set(ls))
            else:
                # Number of distinct values of sum of n_pos items from ls
                distinct_ls = sorted(set(ls))
                possible_sums = set()
                if n_pos == 2:
                    for a in distinct_ls:
                        for b in distinct_ls:
                            possible_sums.add(a + b)
                elif n_pos == 3:
                    for a in distinct_ls:
                        for b in distinct_ls:
                            for c in distinct_ls:
                                possible_sums.add(a + b + c)
                n_classes = len(possible_sums)

            print(f"  {n_pos}-letter {n_items:>7d} {width:>6d} "
                  f"{e_d:>8.4f} {e_d/width:>8.4f} {n_classes:>16d}")
        print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    trigrams = all_trigrams()
    N = len(trigrams)
    assert N == 17576, f"Expected 17576, got {N}"

    print(f"world:// hamming scale aaa-zzz")
    print(f"Trigram count: {N}")
    print(f"Encodings: Binary 15-bit (5+5+5), Ternary 9-trit (3+3+3)")

    for base, wpl, label in [(2, 5, "Binary"), (3, 3, "Ternary")]:
        width = wpl * 3

        print(f"\n{'='*64}")
        print(f"  {label}: {width} {'bits' if base==2 else 'trits'} ({wpl}+{wpl}+{wpl})")
        print(f"{'='*64}")

        # --- 1. Analytical mean ---
        e_d = analytical_mean(base, wpl)
        probs = per_position_disagree_probs(base, wpl)
        print(f"\n  [1] ANALYTICAL MEAN PAIRWISE HAMMING DISTANCE")
        print(f"      Per-digit P(disagree) within one letter:")
        for k, p in enumerate(probs):
            print(f"        digit {k}: {p:.6f}")
        print(f"      Sum per letter: {sum(probs):.6f}")
        print(f"      E[d] for 3-letter trigram: {e_d:.6f}")
        print(f"      E[d] / width: {e_d / width:.6f}")

        # --- 2. Sampled verification ---
        s_mean, s_dists = sampled_mean(trigrams, base, wpl, n_pairs=1000)
        s_counts = Counter(s_dists)
        print(f"\n  [2] SAMPLED VERIFICATION (1000 random pairs)")
        print(f"      Sample mean: {s_mean:.4f}")
        print(f"      Analytical:  {e_d:.4f}")
        print(f"      Difference:  {abs(s_mean - e_d):.4f}")
        print(f"      Distance distribution:")
        for d in range(width + 1):
            c = s_counts.get(d, 0)
            bar = "#" * (c // 5) if c > 0 else ""
            if c > 0:
                print(f"        d={d:2d}: {c:4d}  {bar}")

        # --- 3. Row-sum equivalence classes ---
        ls = letter_row_sums(base, wpl)
        print(f"\n  [3] ROW-SUM EQUIVALENCE CLASSES")
        print(f"      Letter row-sums (single letter, sum to 25 others):")
        distinct_ls = sorted(set(ls))
        ls_counter = Counter(ls)
        for v in distinct_ls:
            letters_with_v = [chr(ord('a') + i) for i, s in enumerate(ls) if s == v]
            print(f"        row-sum={v:3d}: {ls_counter[v]:2d} letters  "
                  f"({','.join(letters_with_v)})")

        sum_counts = compute_row_sum_classes(ls)
        n_classes = len(sum_counts)
        sizes = sorted(sum_counts.values(), reverse=True)
        print(f"\n      Trigram row-sum classes: {n_classes}")
        print(f"      Total trigrams: {sum(sizes)} (check={N})")
        print(f"      Class size distribution:")
        size_dist = Counter(sizes)
        for sz in sorted(size_dist.keys(), reverse=True):
            print(f"        size {sz:>6d}: {size_dist[sz]:3d} classes")

        # Show a few example classes
        sorted_classes = sorted(sum_counts.items(), key=lambda x: x[0])
        print(f"      Example classes (by inner sum):")
        show_indices = [0, len(sorted_classes)//2, -1]
        n2 = 26 ** 2
        for idx in show_indices:
            inner_sum, count = sorted_classes[idx]
            row_sum_val = n2 * inner_sum
            print(f"        inner_sum={inner_sum:4d} → row_sum={row_sum_val:>8d}  "
                  f"  #trigrams={count}")

        # --- 4. Fibonacci-indexed path ---
        fib_idx = fibonacci_indices(N)
        print(f"\n  [4] FIBONACCI-INDEXED TRIGRAM SUBSET")
        print(f"      Fibonacci indices < {N}: {len(fib_idx)} trigrams")
        print(f"      Indices: {fib_idx}")
        fib_trigrams = [trigrams[i] for i in fib_idx]
        print(f"      Trigrams: {fib_trigrams}")

        fwd = path_cost(trigrams, fib_idx, base, wpl)
        rev = path_cost(trigrams, list(reversed(fib_idx)), base, wpl)
        print(f"      Path cost forward:  {fwd}")
        print(f"      Path cost reverse:  {rev}")
        print(f"      Forward == Reverse: {fwd == rev}")

        # Show consecutive distances
        print(f"      Consecutive distances:")
        for k in range(len(fib_idx) - 1):
            a = encode_trigram(trigrams[fib_idx[k]], base, wpl)
            b = encode_trigram(trigrams[fib_idx[k + 1]], base, wpl)
            d = hamming(a, b)
            print(f"        {trigrams[fib_idx[k]]} → {trigrams[fib_idx[k+1]]}  "
                  f"(idx {fib_idx[k]:>5d} → {fib_idx[k+1]:>5d})  d={d}")

    # --- 5. Scaling analysis ---
    scaling_analysis()

    print(f"\n{'='*64}")
    print(f"  DONE")
    print(f"{'='*64}")


if __name__ == "__main__":
    main()
