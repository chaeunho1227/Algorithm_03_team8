"""
데이터셋 일괄 생성 스크립트.

실험 설정:
  N: 200K, 600K, 1.8M
  L: 100
  M: N/10 (10x coverage), N/5 (20x coverage)
  D: 3
  변이율: 1%, 3%, 5%
"""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

from utilis.snp_generator import inject_snps, save_sample, save_snp_log
from utilis.reads_generator import generate_reads, save_reads

BASE = os.path.dirname(__file__)
CHR1 = os.path.join(BASE, 'data/ch1/chr1_chunk.festa')

L    = 100
SEED = 1212

DATASETS = [
    {'n': 200_000,   'label': '200K',  'start': 0},
    {'n': 600_000,   'label': '600K',  'start': 0},
    {'n': 1_800_000, 'label': '1800K', 'start': 0},
]

SNP_RATES = [
    {'label': 'snp_1', 'rate': 0.01},
    {'label': 'snp_3', 'rate': 0.03},
    {'label': 'snp_5', 'rate': 0.05},
]


def load_genome(path: str) -> str:
    sequence = []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence.append(line.strip())
    return ''.join(sequence)


def save_reference(genome: str, path: str, header: str, line_width: int = 80):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'>{header}\n')
        for i in range(0, len(genome), line_width):
            f.write(genome[i:i+line_width] + '\n')
    print(f"Saved: {path} ({len(genome):,} bp)")


def main():
    chr1 = load_genome(CHR1)

    for ds in DATASETS:
        n     = ds['n']
        label = ds['label']
        start = ds['start']
        end   = start + n
        ds_dir   = os.path.join(BASE, f'data/{label}')
        ref_path = os.path.join(ds_dir, f'reference_{label}.festa')

        # reference 추출
        if not os.path.exists(ref_path):
            print(f'\n[{label}] Extracting reference...')
            save_reference(chr1[start:end], ref_path, header=f'reference_{label}')
        else:
            print(f'\n[{label}] Reference already exists, skipping.')

        reference = load_genome(ref_path)

        for snp in SNP_RATES:
            snp_dir     = os.path.join(ds_dir, snp['label'])
            sample_path = os.path.join(snp_dir, 'sample.festa')
            log_path    = os.path.join(snp_dir, 'snp.tsv')

            # sample 생성
            if not os.path.exists(sample_path):
                print(f'  [{snp["label"]}] Generating sample (rate={snp["rate"]})...')
                sample, snp_log = inject_snps(reference, snp_rate=snp['rate'], seed=SEED)
                save_sample(sample, sample_path)
                save_snp_log(snp_log, log_path)
            else:
                print(f'  [{snp["label"]}] Sample already exists, skipping.')
                sample = load_genome(sample_path)

            # reads 생성 (10x, 20x)
            for coverage, m in [(10, n // L * 10), (20, n // L * 20)]:
                reads_path = os.path.join(snp_dir, f'reads_{m}.txt')
                if not os.path.exists(reads_path):
                    print(f'    [reads_{m}] Generating {m:,} reads (L={L}, {coverage}x)...')
                    reads = generate_reads(sample, m=m, l=L, seed=SEED)
                    save_reads(reads, reads_path)
                else:
                    print(f'    [reads_{m}] Already exists, skipping.')

    print('\nAll datasets ready.')


if __name__ == '__main__':
    main()
