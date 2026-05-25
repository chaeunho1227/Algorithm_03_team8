import random
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def reverse_complement(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def generate_reads(genome: str, m: int, l: int, seed: int = None) -> list[tuple[int, str, str]]:
    """
    genome에서 m개 리드를 샘플링한다.

    Returns
    -------
    list of (ref_start, seq, strand)
    """
    rng = random.Random(seed)

    reads = []
    max_start = len(genome) - l

    for _ in range(m):
        pos = rng.randint(0, max_start)
        strand = rng.choice(['+', '-'])
        seq = genome[pos:pos + l]
        if strand == '-':
            seq = reverse_complement(seq)
        reads.append((pos, seq, strand))

    return reads


def save_reads(reads: list[tuple[int, str, str]], out_path: str):
    """reads 서열을 텍스트 파일로 저장 (한 줄 = 한 리드)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True) if os.path.dirname(out_path) else None
    with open(out_path, 'w') as f:
        for _, seq, _ in reads:
            f.write(seq + '\n')
    print(f"Saved: {out_path} ({len(reads):,} reads)")


def save_truth(reads: list[tuple[int, str, str]], out_path: str):
    """정답 매핑 정보를 TSV로 저장 (read별 시작 위치·가닥)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True) if os.path.dirname(out_path) else None
    with open(out_path, 'w') as f:
        f.write('read_name\tref_start\tstrand\n')
        for idx, (pos, _, strand) in enumerate(reads):
            f.write(f'read_{idx:06d}\t{pos}\t{strand}\n')
    print(f"Truth:  {out_path} ({len(reads):,} entries)")