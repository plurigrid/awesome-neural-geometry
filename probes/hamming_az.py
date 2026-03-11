"""
world:// hamming a-z probe — binary (5-bit) and true trits (3-trit)

Explores the Hamming geometry of the 26 lowercase letters (a-z)
under two encodings:

  Binary: 5 bits  (2^5 = 32 ≥ 26)
  Ternary: 3 trits (3^3 = 27 ≥ 26), values {0, 1, 2}

For each encoding:
  1. Codeword representations
  2. Pairwise Hamming distances
  3. Row sums — total Hamming distance from each letter to all others
  4. Subset sums — sums over selected subgroups
  5. Per-position variation analysis

Comparing bit vs trit geometries reveals how the choice of radix
reshapes the metric space: trits pack letters into fewer positions
(3 vs 5) but with a denser alphabet per position (3 vs 2).
"""

import itertools
from collections import Counter


def to_digits(val: int, base: int, width: int) -> tuple[int, ...]:
    """Convert val to a fixed-width tuple of digits in the given base."""
    digits = []
    for _ in range(width):
        digits.append(val % base)
        val //= base
    return tuple(reversed(digits))


def hamming(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    """Hamming distance between two symbol tuples."""
    return sum(x != y for x, y in zip(a, b))


def analyze_encoding(letters, codes, name, base, width):
    """Run full Hamming analysis for a given encoding."""
    n = len(letters)

    print(f"\n{'='*60}")
    print(f"  world:// hamming a-z — {name}")
    print(f"  base={base}  width={width}  capacity={base**width}  used=26")
    print(f"{'='*60}\n")

    # --- Codewords ---
    print(f"Letter  Index  {name}({width}-{'trit' if base==3 else 'bit'})")
    print("-" * 36)
    for ch in letters:
        idx = ord(ch) - ord('a')
        cw = ''.join(str(d) for d in codes[ch])
        print(f"  {ch}       {idx:2d}     {cw}")

    # --- Distance matrix ---
    dist = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            dist[i][j] = hamming(codes[letters[i]], codes[letters[j]])

    print(f"\n--- Pairwise Hamming Distance Matrix ({name}) ---\n")
    header = "    " + " ".join(f"{ch}" for ch in letters)
    print(header)
    print("    " + "--" * n)
    for i, ch in enumerate(letters):
        row = " ".join(f"{dist[i][j]}" for j in range(n))
        print(f" {ch}  {row}")

    # --- Row sums ---
    print(f"\n--- Row Sums ({name}) ---\n")
    row_sums = {}
    for i, ch in enumerate(letters):
        row_sums[ch] = sum(dist[i][j] for j in range(n))
    for ch in letters:
        bar = "#" * (row_sums[ch] // 2)
        print(f"  {ch}: {row_sums[ch]:3d}  {bar}")

    grand = sum(row_sums.values()) // 2
    mean = sum(row_sums.values()) / (n * (n - 1))
    print(f"\n  Grand total pairwise distance: {grand}")
    print(f"  Mean pairwise distance: {mean:.4f}")
    print(f"  Max possible per pair: {width}  "
          f"  Density: {mean/width:.4f}")

    # --- Distance distribution ---
    print(f"\n--- Distance Distribution ({name}) ---\n")
    all_dists = [dist[i][j] for i in range(n) for j in range(i+1, n)]
    dist_counts = Counter(all_dists)
    for d in range(width + 1):
        c = dist_counts.get(d, 0)
        pct = 100 * c / len(all_dists)
        bar = "#" * int(pct / 2)
        print(f"  d={d}: {c:4d} pairs ({pct:5.1f}%)  {bar}")

    # --- Subset sums ---
    print(f"\n--- Subset Sums ({name}) ---\n")
    subsets = {
        "vowels":       ['a', 'e', 'i', 'o', 'u'],
        "first-5":      ['a', 'b', 'c', 'd', 'e'],
        "last-5":       ['v', 'w', 'x', 'y', 'z'],
        "primes":       [chr(ord('a') + p) for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]],
        "powers-of-2":  [chr(ord('a') + p) for p in [1, 2, 4, 8, 16]],
        "fibonacci":    [chr(ord('a') + f) for f in [1, 1, 2, 3, 5, 8, 13, 21]],
    }
    for sname, subset in subsets.items():
        unique = sorted(set(subset))
        total = 0
        count = 0
        for a, b in itertools.combinations(unique, 2):
            ia, ib = ord(a) - ord('a'), ord(b) - ord('a')
            total += dist[ia][ib]
            count += 1
        print(f"  {sname:14s}  letters={','.join(unique):20s}  "
              f"pairs={count:3d}  sum={total:4d}  "
              f"mean={total/max(count,1):.3f}")

    # --- Per-position analysis ---
    print(f"\n--- Per-Position Analysis ({name}) ---\n")
    for pos in range(width):
        vals = [codes[ch][pos] for ch in letters]
        val_counts = Counter(vals)
        # pairs differing at this position = total_pairs - pairs_same
        pairs_same = sum(c * (c - 1) // 2 for c in val_counts.values())
        pairs_diff = n * (n - 1) // 2 - pairs_same
        dist_str = " ".join(f"{s}:{val_counts.get(s,0):2d}" for s in range(base))
        print(f"  Pos {pos} ({base}^{width-1-pos}={base**(width-1-pos):2d}):  "
              f"[{dist_str}]  pairs differing={pairs_diff:3d}")

    total_pos_diff = sum(
        codes[letters[i]][p] != codes[letters[j]][p]
        for i in range(n) for j in range(i+1, n)
        for p in range(width)
    )
    print(f"\n  Total position-disagreements: {total_pos_diff}")
    print(f"  (Equals grand total pairwise Hamming: {grand})")

    return dist, row_sums, grand, mean


def compare(binary_stats, ternary_stats):
    """Compare binary vs ternary Hamming geometries."""
    _, b_sums, b_grand, b_mean = binary_stats
    _, t_sums, t_grand, t_mean = ternary_stats
    letters = sorted(b_sums.keys())

    print(f"\n{'='*60}")
    print(f"  COMPARISON: Binary (5-bit) vs True Trits (3-trit)")
    print(f"{'='*60}\n")

    print(f"  {'Metric':<35s} {'Binary':>8s} {'Ternary':>8s} {'Ratio':>8s}")
    print(f"  {'-'*35} {'-'*8} {'-'*8} {'-'*8}")
    print(f"  {'Positions per codeword':<35s} {'5':>8s} {'3':>8s} {'0.60':>8s}")
    print(f"  {'Symbols per position':<35s} {'2':>8s} {'3':>8s} {'1.50':>8s}")
    print(f"  {'Codespace capacity':<35s} {'32':>8s} {'27':>8s} {'0.84':>8s}")
    print(f"  {'Unused codewords':<35s} {'6':>8s} {'1':>8s} {'0.17':>8s}")
    print(f"  {'Grand total pairwise distance':<35s} {b_grand:>8d} {t_grand:>8d} "
          f"{t_grand/b_grand:>8.4f}")
    print(f"  {'Mean pairwise distance':<35s} {b_mean:>8.4f} {t_mean:>8.4f} "
          f"{t_mean/b_mean:>8.4f}")
    print(f"  {'Mean / max (density)':<35s} {b_mean/5:>8.4f} {t_mean/3:>8.4f} "
          f"{(t_mean/3)/(b_mean/5):>8.4f}")

    # Per-letter comparison
    print(f"\n--- Per-Letter Row Sums: Binary vs Ternary ---\n")
    print(f"  Letter  Binary  Ternary  Diff")
    print(f"  {'-'*6}  {'-'*6}  {'-'*7}  {'-'*4}")
    for ch in letters:
        diff = t_sums[ch] - b_sums[ch]
        sign = "+" if diff > 0 else " " if diff == 0 else ""
        print(f"    {ch}      {b_sums[ch]:4d}    {t_sums[ch]:5d}   {sign}{diff}")

    # Rank correlation
    b_ranked = sorted(letters, key=lambda c: b_sums[c])
    t_ranked = sorted(letters, key=lambda c: t_sums[c])
    b_rank = {ch: i for i, ch in enumerate(b_ranked)}
    t_rank = {ch: i for i, ch in enumerate(t_ranked)}
    n = len(letters)
    d_sq = sum((b_rank[ch] - t_rank[ch]) ** 2 for ch in letters)
    spearman = 1 - 6 * d_sq / (n * (n * n - 1))
    print(f"\n  Spearman rank correlation of row sums: {spearman:.4f}")

    # Which letters change rank the most?
    rank_shifts = [(ch, abs(b_rank[ch] - t_rank[ch])) for ch in letters]
    rank_shifts.sort(key=lambda x: -x[1])
    print(f"\n  Largest rank shifts (binary→ternary):")
    for ch, shift in rank_shifts[:5]:
        print(f"    {ch}: rank {b_rank[ch]:2d} → {t_rank[ch]:2d}  (shift={shift})")


def main():
    letters = [chr(ord('a') + i) for i in range(26)]

    # Binary: 5 bits
    bits = {ch: to_digits(ord(ch) - ord('a'), 2, 5) for ch in letters}
    binary_stats = analyze_encoding(letters, bits, "Binary", 2, 5)

    # True Trits: 3 trits
    trits = {ch: to_digits(ord(ch) - ord('a'), 3, 3) for ch in letters}
    ternary_stats = analyze_encoding(letters, trits, "True Trits", 3, 3)

    # Comparison
    compare(binary_stats, ternary_stats)


if __name__ == "__main__":
    main()
