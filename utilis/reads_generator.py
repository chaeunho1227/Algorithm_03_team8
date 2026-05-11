import random
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from utilis.extract_region import load_genome_raw


def generate_reads(genome: str, m: int, l: int, seed: int = None) -> list[tuple[int, str]]:
    if seed is not None:
        random.seed(seed)

    reads = []
    max_start = len(genome) - l
    positions = random.choices(range(max_start), k=m)

    for pos in positions:
        reads.append((pos, genome[pos:pos+l]))

    return reads


def save_reads(reads: list[tuple[int, str]], out_path: str):
    """reads를 텍스트 파일로 저장. 각 줄에 read 서열."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True) if os.path.dirname(out_path) else None
    with open(out_path, 'w') as f:
        for _, read in reads:
            f.write(read + '\n')
    print(f"Saved: {out_path} ({len(reads):,} reads)")