from __future__ import annotations
import random
import os


def generate_reads(genome: str, m: int, l: int, seed: int = None) -> list[tuple[int, str]]:
    rng = random.Random(seed)
    reads = []
    max_start = len(genome) - l

    for _ in range(m):
        pos = rng.randint(0, max_start)
        reads.append((pos, genome[pos:pos + l]))

    return reads


def save_reads(reads: list[tuple[int, str]], out_path: str):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        for _, seq in reads:
            f.write(seq + '\n')
    print(f"Saved: {out_path} ({len(reads):,} reads)")


def save_truth(reads: list[tuple[int, str]], out_path: str):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write('read_name\tref_start\n')
        for idx, (pos, _) in enumerate(reads):
            f.write(f'read_{idx:06d}\t{pos}\n')
    print(f"Truth:  {out_path} ({len(reads):,} entries)")
