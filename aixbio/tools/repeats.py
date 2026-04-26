from __future__ import annotations


def find_direct_repeats(dna: str, min_len: int = 20) -> list[tuple[int, int, int]]:
    """Find non-overlapping direct repeats of length >= min_len.

    Returns a list of (pos1, pos2, length) tuples.  Long perfect direct repeats
    (>= 20 bp) can trigger RecA-independent deletion in E. coli during plasmid
    replication, making the construct unstable.

    Algorithm: index all k-mers at min_len, then extend each pair to its
    maximal non-overlapping match.  O(n * min_len) — suitable up to ~10 kb.
    """
    dna = dna.upper()
    n = len(dna)
    if n < min_len * 2:
        return []

    kmer_positions: dict[str, list[int]] = {}
    for i in range(n - min_len + 1):
        kmer = dna[i : i + min_len]
        kmer_positions.setdefault(kmer, []).append(i)

    results: list[tuple[int, int, int]] = []
    seen: set[tuple[int, int]] = set()

    for positions in kmer_positions.values():
        if len(positions) < 2:
            continue
        for a_idx in range(len(positions)):
            for b_idx in range(a_idx + 1, len(positions)):
                i, j = positions[a_idx], positions[b_idx]
                if (i, j) in seen:
                    continue
                length = min_len
                while j + length < n and dna[i + length] == dna[j + length]:
                    length += 1
                if j >= i + length:
                    seen.add((i, j))
                    results.append((i, j, length))

    results.sort(key=lambda t: -t[2])
    return results