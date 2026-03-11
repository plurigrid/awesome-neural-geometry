"""
world:// hamming aa-zz probe — bigram Hamming geometry (676 bigrams)

Explores the Hamming geometry of all 676 two-letter bigrams (aa..zz)
under two encodings:

  Binary:  10 bits (5 per letter, concatenated)
  Ternary:  6 trits (3 per letter, concatenated)

Computes:
  1. Distance-profile equivalence classes (bigrams sharing identical
     sorted distance rows) — count and sizes only
  2. Grand total pairwise Hamming distance, mean, density (mean/max)
  3. Fibonacci-indexed bigram subset path costs (forward & reverse),
     palindrome-invariance check
  4. Number of "paren colors" needed (= number of profile equiv classes)

The 676x676 distance matrix is computed but NOT printed.
"""

import itertools
from collections import Counter, defaultdict


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


def make_bigrams():
    """Return the 676 bigrams aa..zz in lexicographic order."""
    return [chr(ord('a') + i) + chr(ord('a') + j)
            for i in range(26) for j in range(26)]


def encode_bigram(bg: str, base: int, letter_width: int) -> tuple[int, ...]:
    """Encode a bigram by concatenating per-letter codes."""
    c0 = to_digits(ord(bg[0]) - ord('a'), base, letter_width)
    c1 = to_digits(ord(bg[1]) - ord('a'), base, letter_width)
    return c0 + c1


def fibonacci_indices(limit: int) -> list[int]:
    """Return Fibonacci numbers (starting 1,2,3,5,...) that are < limit."""
    fibs = []
    a, b = 1, 2
    while a < limit:
        fibs.append(a)
        a, b = b, a + b
    return fibs


def analyze_bigram_encoding(bigrams, codes, name, base, letter_width):
    """Run bigram Hamming analysis for one encoding."""
    n = len(bigrams)
    width = letter_width * 2  # total code width per bigram

    print(f"\n{'='*64}")
    print(f"  world:// hamming aa-zz — {name}")
    print(f"  base={base}  letter_width={letter_width}  "
          f"bigram_width={width}  bigrams={n}")
    print(f"{'='*64}")

    # --- Compute distance rows (not the full matrix, just each row sorted) ---
    # We need: row sums, sorted distance profiles, grand total, fib path
    # For efficiency, compute pairwise distances on the fly.

    # Build code list for fast access
    code_list = [codes[bg] for bg in bigrams]

    # Compute row sums and sorted profiles
    print(f"\n--- Computing {n}x{n} pairwise Hamming distances ---")
    row_sums = [0] * n
    grand_total = 0
    dist_counter = Counter()
    sorted_profiles = []

    for i in range(n):
        row = []
        for j in range(n):
            d = hamming(code_list[i], code_list[j])
            row.append(d)
            if j > i:
                grand_total += d
                dist_counter[d] += 1
            row_sums[i] += d
        sorted_profiles.append(tuple(sorted(row)))

    # --- 1. Distance-profile equivalence classes ---
    print(f"\n--- 1. Distance-Profile Equivalence Classes ({name}) ---\n")
    profile_classes = defaultdict(list)
    for i, prof in enumerate(sorted_profiles):
        profile_classes[prof].append(bigrams[i])

    class_sizes = sorted([len(v) for v in profile_classes.values()], reverse=True)
    num_classes = len(profile_classes)
    print(f"  Number of equivalence classes: {num_classes}")
    size_counts = Counter(class_sizes)
    print(f"  Class size distribution:")
    for sz in sorted(size_counts.keys(), reverse=True):
        print(f"    size {sz:3d}: {size_counts[sz]:3d} classes")
    if len(class_sizes) <= 30:
        print(f"  All class sizes: {class_sizes}")

    # --- 2. Grand total, mean, density ---
    print(f"\n--- 2. Grand Total Pairwise Distance ({name}) ---\n")
    num_pairs = n * (n - 1) // 2
    mean_dist = grand_total / num_pairs
    max_dist = width
    density = mean_dist / max_dist
    print(f"  Grand total pairwise distance: {grand_total}")
    print(f"  Number of pairs:               {num_pairs}")
    print(f"  Mean pairwise distance:        {mean_dist:.6f}")
    print(f"  Max possible per pair:         {max_dist}")
    print(f"  Density (mean/max):            {density:.6f}")

    # Distance distribution
    print(f"\n  Distance distribution:")
    for d in range(max_dist + 1):
        c = dist_counter.get(d, 0)
        pct = 100 * c / num_pairs
        bar = "#" * int(pct / 2)
        print(f"    d={d:2d}: {c:6d} pairs ({pct:5.1f}%)  {bar}")

    # --- 3. Fibonacci-indexed bigram subset ---
    print(f"\n--- 3. Fibonacci-Indexed Bigram Subset ({name}) ---\n")
    fib_idx = fibonacci_indices(n)
    print(f"  Fibonacci indices (<{n}): {fib_idx}")
    fib_bigrams = [bigrams[i] for i in fib_idx]
    print(f"  Bigrams: {' '.join(fib_bigrams)}")

    # Forward path cost
    fwd_cost = 0
    for k in range(len(fib_bigrams) - 1):
        d = hamming(codes[fib_bigrams[k]], codes[fib_bigrams[k + 1]])
        fwd_cost += d
    print(f"  Forward path cost:  {fwd_cost}")

    # Reverse path cost
    rev_bigrams = list(reversed(fib_bigrams))
    rev_cost = 0
    for k in range(len(rev_bigrams) - 1):
        d = hamming(codes[rev_bigrams[k]], codes[rev_bigrams[k + 1]])
        rev_cost += d
    print(f"  Reverse path cost:  {rev_cost}")

    palindrome_inv = fwd_cost == rev_cost
    print(f"  Palindrome-invariant? {palindrome_inv}")

    # Show step-by-step
    print(f"\n  Forward path steps:")
    for k in range(len(fib_bigrams) - 1):
        a, b = fib_bigrams[k], fib_bigrams[k + 1]
        d = hamming(codes[a], codes[b])
        print(f"    {a} -> {b} : d={d}")

    # --- 4. Paren colors needed ---
    print(f"\n--- 4. Paren Colors Needed ({name}) ---\n")
    print(f"  Colors needed = number of profile equiv classes = {num_classes}")

    return num_classes, grand_total, mean_dist, density


def main():
    bigrams = make_bigrams()
    assert len(bigrams) == 676

    # Binary: 5 bits per letter -> 10 bits per bigram
    bin_codes = {bg: encode_bigram(bg, 2, 5) for bg in bigrams}
    bin_stats = analyze_bigram_encoding(bigrams, bin_codes, "Binary (10-bit)", 2, 5)

    # Ternary: 3 trits per letter -> 6 trits per bigram
    tri_codes = {bg: encode_bigram(bg, 3, 3) for bg in bigrams}
    tri_stats = analyze_bigram_encoding(bigrams, tri_codes, "Ternary (6-trit)", 3, 3)

    # --- Comparison ---
    print(f"\n{'='*64}")
    print(f"  COMPARISON: Binary (10-bit) vs Ternary (6-trit)")
    print(f"{'='*64}\n")
    labels = ["Profile equiv classes", "Grand total pairwise",
              "Mean pairwise dist", "Density (mean/max)"]
    for label, bv, tv in zip(labels, bin_stats, tri_stats):
        if isinstance(bv, float):
            print(f"  {label:<28s}  Binary: {bv:12.6f}  Ternary: {tv:12.6f}")
        else:
            print(f"  {label:<28s}  Binary: {bv:12}  Ternary: {tv:12}")


if __name__ == "__main__":
    main()
