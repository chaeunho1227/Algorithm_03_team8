"""
랜덤 reference / sample / reads / truth 데이터셋 생성 스크립트 (제출용)

3개 genome 크기 × 3개 SNP rate × 2개 coverage 조합으로 데이터셋을 생성한다.
이미 존재하는 파일은 건너뛴다.

실행:
    python3 generate_random_datasets.py

출력 디렉토리 구조 (DATA_DIR 기준):
    {DATA_DIR}/200K/reference_200K.festa
    {DATA_DIR}/200K/snp_1/sample.festa
    {DATA_DIR}/200K/snp_1/reads_20000.txt
    {DATA_DIR}/200K/snp_1/truth_20000.tsv
    ...
"""

import os
import random


# ==============================================================
#  경로 / 파라미터 설정 (제출 환경에 맞게 수정)
# ==============================================================

# 생성된 데이터셋을 저장할 루트 디렉토리
# run_experiments.py 의 DATA_DIR 과 동일하게 맞춘다.
DATA_DIR = "data"

# read 길이 (bp)
L = 100

# 난수 시드 (재현성)
SEED = 1212

# 테스트할 genome 크기
DATASETS = [
    {"n": 200_000,   "label": "200K"},
    {"n": 600_000,   "label": "600K"},
    {"n": 1_800_000, "label": "1800K"},
]

# 테스트할 SNP rate
SNP_RATES = [
    {"label": "snp_1", "rate": 0.01},
    {"label": "snp_3", "rate": 0.03},
    {"label": "snp_5", "rate": 0.05},
]

# 테스트할 coverage
COVERAGES = [10, 20]


# ==============================================================
#  SNP 주입 / sample 저장
# ==============================================================

def inject_snps(genome, snp_rate, seed=None):
    if seed is not None:
        random.seed(seed)

    genome = list(genome)
    bases = {"A", "C", "G", "T"}

    n_snps = int(len(genome) * snp_rate)
    positions = random.sample(range(len(genome)), n_snps)

    for pos in positions:
        original = genome[pos]
        genome[pos] = random.choice(list(bases - {original}))

    return "".join(genome)


def save_sample(genome, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write(">sample\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + "\n")
    print(f"Saved: {out_path} ({len(genome):,} bp)")


# ==============================================================
#  reads / truth 생성
# ==============================================================

def generate_reads(genome, m, l, seed=None):
    rng = random.Random(seed)
    reads = []
    max_start = len(genome) - l

    for _ in range(m):
        pos = rng.randint(0, max_start)
        reads.append((pos, genome[pos:pos + l]))

    return reads


def save_reads(reads, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        for _, seq in reads:
            f.write(seq + "\n")
    print(f"Saved: {out_path} ({len(reads):,} reads)")


def save_truth(reads, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("read_name\tref_start\n")
        for idx, (pos, _) in enumerate(reads):
            f.write(f"read_{idx:06d}\t{pos}\n")
    print(f"Truth:  {out_path} ({len(reads):,} entries)")


# ==============================================================
#  reference 로드 / 저장
# ==============================================================

def load_genome(path):
    sequence = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip())
    return "".join(sequence)


def save_reference(genome, path, header):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + "\n")
    print(f"Saved: {path} ({len(genome):,} bp)")


# ==============================================================
#  main
# ==============================================================

def main():
    random.seed(SEED)

    for ds in DATASETS:
        n, label = ds["n"], ds["label"]
        ds_dir   = os.path.join(DATA_DIR, label)
        ref_path = os.path.join(ds_dir, f"reference_{label}.festa")

        if not os.path.exists(ref_path):
            print(f"\n[{label}] Generating random reference ({n:,} bp)...")
            genome = "".join(random.choices("ACGT", k=n))
            save_reference(genome, ref_path, header=f"reference_random_{label}")
        else:
            print(f"\n[{label}] Reference already exists, skipping.")

        reference = load_genome(ref_path)

        for snp in SNP_RATES:
            snp_dir     = os.path.join(ds_dir, snp["label"])
            sample_path = os.path.join(snp_dir, "sample.festa")

            if not os.path.exists(sample_path):
                print(f"  [{snp['label']}] Generating sample (rate={snp['rate']})...")
                sample = inject_snps(reference, snp_rate=snp["rate"], seed=SEED)
                save_sample(sample, sample_path)
            else:
                print(f"  [{snp['label']}] Sample already exists, skipping.")
                sample = load_genome(sample_path)

            for coverage in COVERAGES:
                m          = n // L * coverage
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                truth_path = os.path.join(snp_dir, f"truth_{m}.tsv")
                if not os.path.exists(reads_path):
                    print(f"    [reads_{m}] Generating {m:,} reads ({coverage}x)...")
                    reads = generate_reads(sample, m=m, l=L, seed=SEED)
                    save_reads(reads, reads_path)
                    save_truth(reads, truth_path)
                else:
                    print(f"    [reads_{m}] Already exists, skipping.")

    print("\nAll random datasets ready.")


if __name__ == "__main__":
    main()
