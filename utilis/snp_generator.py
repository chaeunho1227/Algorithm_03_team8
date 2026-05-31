from __future__ import annotations
import random
import os


def inject_snps(genome: str, snp_rate: float, seed: int = None) -> str:
    if seed is not None:
        random.seed(seed)

    genome = list(genome)
    bases = {'A', 'C', 'G', 'T'}

    n_snps = int(len(genome) * snp_rate)
    positions = random.sample(range(len(genome)), n_snps)

    for pos in positions:
        original = genome[pos]
        genome[pos] = random.choice(list(bases - {original}))

    return ''.join(genome)


def save_sample(genome: str, out_path: str):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write('>sample\n')
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + '\n')
    print(f"Saved: {out_path} ({len(genome):,} bp)")
