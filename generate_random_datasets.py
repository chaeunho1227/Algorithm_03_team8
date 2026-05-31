import os
import sys
import random
sys.path.insert(0, os.path.dirname(__file__))

from utilis.snp_generator import inject_snps, save_sample
from utilis.reads_generator import generate_reads, save_reads, save_truth

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'random')

L    = 100
SEED = 1212

DATASETS = [
    {'n': 200_000,   'label': '200K'},
    {'n': 600_000,   'label': '600K'},
    {'n': 1_800_000, 'label': '1800K'},
]

SNP_RATES = [
    {'label': 'snp_1', 'rate': 0.01},
    {'label': 'snp_3', 'rate': 0.03},
    {'label': 'snp_5', 'rate': 0.05},
]


def load_genome(path):
    sequence = []
    with open(path) as f:
        for line in f:
            if not line.startswith('>'):
                sequence.append(line.strip())
    return ''.join(sequence)


def save_reference(genome, path, header):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'>{header}\n')
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + '\n')
    print(f"Saved: {path} ({len(genome):,} bp)")


def main():
    random.seed(SEED)

    for ds in DATASETS:
        n, label = ds['n'], ds['label']
        ds_dir   = os.path.join(DATA_DIR, label)
        ref_path = os.path.join(ds_dir, f'reference_{label}.festa')

        if not os.path.exists(ref_path):
            print(f'\n[random/{label}] Generating random reference ({n:,} bp)...')
            genome = ''.join(random.choices('ACGT', k=n))
            save_reference(genome, ref_path, header=f'reference_random_{label}')
        else:
            print(f'\n[random/{label}] Reference already exists, skipping.')

        reference = load_genome(ref_path)

        for snp in SNP_RATES:
            snp_dir     = os.path.join(ds_dir, snp['label'])
            sample_path = os.path.join(snp_dir, 'sample.festa')

            if not os.path.exists(sample_path):
                print(f'  [{snp["label"]}] Generating sample (rate={snp["rate"]})...')
                sample = inject_snps(reference, snp_rate=snp['rate'], seed=SEED)
                save_sample(sample, sample_path)
            else:
                print(f'  [{snp["label"]}] Sample already exists, skipping.')
                sample = load_genome(sample_path)

            for coverage, m in [(10, n // L * 10), (20, n // L * 20)]:
                reads_path = os.path.join(snp_dir, f'reads_{m}.txt')
                truth_path = os.path.join(snp_dir, f'truth_{m}.tsv')
                if not os.path.exists(reads_path):
                    print(f'    [reads_{m}] Generating {m:,} reads ({coverage}x)...')
                    reads = generate_reads(sample, m=m, l=L, seed=SEED)
                    save_reads(reads, reads_path)
                    save_truth(reads, truth_path)
                else:
                    print(f'    [reads_{m}] Already exists, skipping.')

    print('\nAll random datasets ready.')


if __name__ == '__main__':
    main()
