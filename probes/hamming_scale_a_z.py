"""
world:// hamming scale a-z — distance matrices, equivalence classes, coloring

Computes for all 26 letters (a-z) under two encodings:
  Binary:  5 bits  (indices 0-25 in base 2, width 5)
  Ternary: 3 trits (indices 0-25 in base 3, width 3)

Reports:
  1. Codewords for each letter in both encodings
  2. Full 26x26 Hamming distance matrices (computed, not printed)
  3. Distance-profile equivalence classes (letters sharing identical
     sorted distance rows)
  4. Number of distinct paren-coloring classes needed
  5. Summary statistics: grand total, mean, density
"""

from collections import defaultdict


def to_digits(val: int, base: int, width: int) -> tuple[int, ...]:
    """Convert val to a fixed-width tuple of digits in the given base."""
    digits = []
    for _ in range(width):
        digits.append(val % base)
        val //= base
    return tuple(reversed(digits))


def hamming(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x != y for x, y in zip(a, b))


def build_distance_matrix(letters, codes):
    """Return 26x26 Hamming distance matrix."""
    n = len(letters)
    dist = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = hamming(codes[letters[i]], codes[letters[j]])
            dist[i][j] = d
            dist[j][i] = d
    return dist


def equivalence_classes(letters, dist):
    """Group letters by their sorted distance profile (row of distances to all others)."""
    profiles = defaultdict(list)
    for i, ch in enumerate(letters):
        # The sorted row (excluding self-distance 0) is the profile signature
        profile = tuple(sorted(dist[i][j] for j in range(len(letters)) if j != i))
        profiles[profile].append(ch)
    return profiles


def summary_stats(letters, dist, width):
    """Compute grand total, mean pairwise distance, and density."""
    n = len(letters)
    total = sum(dist[i][j] for i in range(n) for j in range(i + 1, n))
    pairs = n * (n - 1) // 2
    mean = total / pairs
    density = mean / width
    return total, mean, density


def report(name, base, width, letters, codes, dist):
    """Print compact report for one encoding."""
    n = len(letters)

    print(f"\n{'='*60}")
    print(f"  {name}: base={base}  width={width}  "
          f"capacity={base**width}  used={n}")
    print(f"{'='*60}")

    # Codewords
    print(f"\n  Codewords ({name}):")
    rows = []
    for ch in letters:
        cw = ''.join(str(d) for d in codes[ch])
        rows.append(f"{ch}={cw}")
    # Print 13 per line for compactness
    for start in range(0, len(rows), 13):
        print(f"    {' '.join(rows[start:start+13])}")

    # Summary stats
    total, mean, density = summary_stats(letters, dist, width)
    print(f"\n  Grand total pairwise distance : {total}")
    print(f"  Mean pairwise distance       : {mean:.4f}")
    print(f"  Max possible per pair         : {width}")
    print(f"  Density (mean/max)            : {density:.4f}")

    # Equivalence classes
    classes = equivalence_classes(letters, dist)
    print(f"\n  Distance-profile equivalence classes: {len(classes)}")
    for idx, (profile, members) in enumerate(
        sorted(classes.items(), key=lambda kv: (-len(kv[1]), kv[1][0]))
    ):
        profile_hist = {}
        for d in profile:
            profile_hist[d] = profile_hist.get(d, 0) + 1
        hist_str = " ".join(f"d{d}:{c}" for d, c in sorted(profile_hist.items()))
        print(f"    Class {idx}: [{','.join(members)}]  "
              f"(size {len(members)})  dist-hist: {hist_str}")

    # Coloring classes needed
    print(f"\n  Distinct paren-coloring classes needed: {len(classes)}")

    return total, mean, density, classes


def main():
    letters = [chr(ord('a') + i) for i in range(26)]

    # --- Binary encoding: 5 bits ---
    bin_codes = {ch: to_digits(ord(ch) - ord('a'), 2, 5) for ch in letters}
    bin_dist = build_distance_matrix(letters, bin_codes)
    b_total, b_mean, b_dens, b_classes = report(
        "Binary (5-bit)", 2, 5, letters, bin_codes, bin_dist
    )

    # --- Ternary encoding: 3 trits ---
    tri_codes = {ch: to_digits(ord(ch) - ord('a'), 3, 3) for ch in letters}
    tri_dist = build_distance_matrix(letters, tri_codes)
    t_total, t_mean, t_dens, t_classes = report(
        "Ternary (3-trit)", 3, 3, letters, tri_codes, tri_dist
    )

    # --- Cross-encoding comparison ---
    print(f"\n{'='*60}")
    print(f"  Cross-encoding comparison")
    print(f"{'='*60}")
    print(f"\n  {'Metric':<35s} {'Binary':>10s} {'Ternary':>10s}")
    print(f"  {'-'*35} {'-'*10} {'-'*10}")
    print(f"  {'Grand total distance':<35s} {b_total:>10d} {t_total:>10d}")
    print(f"  {'Mean pairwise distance':<35s} {b_mean:>10.4f} {t_mean:>10.4f}")
    print(f"  {'Density (mean/max)':<35s} {b_dens:>10.4f} {t_dens:>10.4f}")
    print(f"  {'Equivalence classes':<35s} {len(b_classes):>10d} {len(t_classes):>10d}")
    print(f"  {'Coloring classes needed':<35s} {len(b_classes):>10d} {len(t_classes):>10d}")
    print()


if __name__ == "__main__":
    main()
