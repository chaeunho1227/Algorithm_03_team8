import random
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def inject_snps(genome: str, snp_rate: float, seed: int = None) -> tuple[str, list[tuple[int, str, str]]]:
    if seed is not None:
        random.seed(seed)

    genome = list(genome)
    snp_log = []
    bases = {'A', 'C', 'G', 'T'}

    n_snps = int(len(genome) * snp_rate)
    positions = random.sample(range(len(genome)), n_snps)

    for pos in positions:
        original = genome[pos]
        mutated = random.choice(list(bases - {original}))
        genome[pos] = mutated
        snp_log.append((pos, original, mutated))

    snp_log.sort(key=lambda x: x[0])
    return ''.join(genome), snp_log


def save_sample(genome: str, out_path: str, header: str = 'sample', line_width: int = 80):
    os.makedirs(os.path.dirname(out_path), exist_ok=True) if os.path.dirname(out_path) else None
    with open(out_path, 'w') as f:
        f.write(f'>{header}\n')
        for i in range(0, len(genome), line_width):
            f.write(genome[i:i+line_width] + '\n')
    print(f"Saved: {out_path} ({len(genome):,} bp)")


def save_snp_log(snp_log: list, out_path: str):
    os.makedirs(os.path.dirname(out_path), exist_ok=True) if os.path.dirname(out_path) else None
    with open(out_path, 'w') as f:
        f.write('position\toriginal\tmutated\n')
        for pos, orig, mut in snp_log:
            f.write(f'{pos}\t{orig}\t{mut}\n')
    print(f"SNP log: {out_path} ({len(snp_log):,} SNPs)")