"""
world:// hamming a-z probe

Explores the Hamming geometry of the 26 lowercase letters (a-z)
by treating each letter as its ASCII binary codeword and computing:

  1. Binary representations (5-bit, since a-z spans 0..25)
  2. Pairwise Hamming distances
  3. Row sums — total Hamming distance from each letter to all others
  4. Column/any sums — sums over arbitrary subsets

This connects to neural geometry through the lens of discrete
metric spaces, error-correcting codes, and binary neural codes.
"""

import itertools


def letter_to_bits(ch: str, width: int = 5) -> tuple[int, ...]:
    """Map a lowercase letter to its binary tuple (0-indexed: a=0, b=1, ..., z=25)."""
    val = ord(ch) - ord('a')
    return tuple((val >> (width - 1 - i)) & 1 for i in range(width))


def hamming(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    """Hamming distance between two binary tuples."""
    return sum(x != y for x, y in zip(a, b))


def main():
    letters = [chr(ord('a') + i) for i in range(26)]
    bits = {ch: letter_to_bits(ch) for ch in letters}

    # --- Binary representations ---
    print("=== world:// hamming a-z ===\n")
    print("Letter  Index  Binary(5-bit)")
    print("-" * 32)
    for ch in letters:
        idx = ord(ch) - ord('a')
        print(f"  {ch}       {idx:2d}     {''.join(str(b) for b in bits[ch])}")

    # --- Pairwise Hamming distance matrix ---
    n = len(letters)
    dist = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            dist[i][j] = hamming(bits[letters[i]], bits[letters[j]])

    print("\n=== Pairwise Hamming Distance Matrix (26x26) ===\n")
    header = "    " + " ".join(f"{ch}" for ch in letters)
    print(header)
    print("    " + "--" * n)
    for i, ch in enumerate(letters):
        row = " ".join(f"{dist[i][j]}" for j in range(n))
        print(f" {ch}  {row}")

    # --- Row sums: total Hamming distance from each letter to all others ---
    print("\n=== Row Sums (total distance from each letter to all others) ===\n")
    row_sums = {}
    for i, ch in enumerate(letters):
        s = sum(dist[i][j] for j in range(n))
        row_sums[ch] = s
    for ch in letters:
        bar = "#" * (row_sums[ch] // 2)
        print(f"  {ch}: {row_sums[ch]:3d}  {bar}")

    print(f"\n  Grand total (sum of all pairwise distances): {sum(row_sums.values()) // 2}")
    print(f"  Mean pairwise distance: {sum(row_sums.values()) / (n * (n - 1)):.4f}")

    # --- Sums over subsets ("sums only of any of them") ---
    print("\n=== Subset Hamming Sums (sums of distances within selected subsets) ===\n")

    interesting_subsets = {
        "vowels":       ['a', 'e', 'i', 'o', 'u'],
        "first-5":      ['a', 'b', 'c', 'd', 'e'],
        "last-5":       ['v', 'w', 'x', 'y', 'z'],
        "primes":       [chr(ord('a') + p) for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]],
        "powers-of-2":  [chr(ord('a') + p) for p in [1, 2, 4, 8, 16]],
        "fibonacci":    [chr(ord('a') + f) for f in [1, 1, 2, 3, 5, 8, 13, 21]],
    }

    for name, subset in interesting_subsets.items():
        unique = sorted(set(subset))
        total = 0
        count = 0
        for a, b in itertools.combinations(unique, 2):
            ia, ib = ord(a) - ord('a'), ord(b) - ord('a')
            total += dist[ia][ib]
            count += 1
        print(f"  {name:14s}  letters={','.join(unique):20s}  "
              f"pairs={count:3d}  sum={total:4d}  "
              f"mean={total/max(count,1):.3f}")

    # --- Per-bit statistics ---
    print("\n=== Per-Bit Analysis (which bit positions carry the most variation) ===\n")
    for bit_pos in range(5):
        ones = sum(bits[ch][bit_pos] for ch in letters)
        zeros = 26 - ones
        contribution = ones * zeros  # number of pairs differing at this bit
        print(f"  Bit {bit_pos} (2^{4-bit_pos}={1<<(4-bit_pos):2d}):  "
              f"zeros={zeros:2d}  ones={ones:2d}  "
              f"pairs differing={contribution:3d}")

    total_pairs_differing = sum(
        bits[letters[i]][bp] != bits[letters[j]][bp]
        for i in range(n) for j in range(i+1, n)
        for bp in range(5)
    )
    print(f"\n  Total bit-disagreements across all pairs: {total_pairs_differing}")
    print(f"  (This equals the grand total of pairwise Hamming distances: "
          f"{sum(row_sums.values()) // 2})")


if __name__ == "__main__":
    main()
