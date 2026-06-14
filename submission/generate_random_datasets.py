"""랜덤 서열 기반 데이터셋 생성.

실제 chr1 reference 대신 ACGT를 무작위로 생성한 reference를 사용한다.
랜덤 서열은 반복 구간이 거의 없어 k-mer가 고유해지므로, 매핑 알고리즘
자체의 성능(정렬률·복원정확도)을 반복서열 노이즈 없이 깨끗하게 측정할 수 있다.

각 genome 크기 × SNP rate × coverage 조합마다 다음을 생성한다:
- reference_{label}.festa : 원본 reference (랜덤 ACGT)
- sample.festa           : reference에 SNP를 주입한 정답 게놈
- reads_{m}.txt          : sample에서 무작위로 떼어낸 read들
- truth_{m}.tsv          : 각 read의 정답 시작 위치
"""
import os
import random


DATA_DIR = "data"
L = 100        # read 한 개의 길이 (bp)
SEED = 1212    # 재현성 확보용 고정 시드

# 테스트할 genome 크기 (200K / 600K / 1800K)
DATASETS = [
    {"n": 200_000, "label": "200K"},
    {"n": 600_000, "label": "600K"},
    {"n": 1_800_000, "label": "1800K"},
]

# 테스트할 SNP rate (1% / 3% / 5%)
SNP_RATES = [
    {"label": "snp_1", "rate": 0.01},
    {"label": "snp_3", "rate": 0.03},
    {"label": "snp_5", "rate": 0.05},
]

# 테스트할 coverage (10x / 20x)
COVERAGES = [10, 20]


def inject_snps(genome, snp_rate, seed=None):
    """reference에 snp_rate 비율만큼 SNP(치환 변이)를 주입한 sample 게놈을 생성.

    Returns: SNP가 적용된 게놈 문자열 (reference와 동일 길이).
    """
    if seed is not None:
        random.seed(seed)

    genome = list(genome)
    bases = {"A", "C", "G", "T"}

    # 변이시킬 위치를 중복 없이(sample) 무작위 선택한다.
    n_snps = int(len(genome) * snp_rate)
    positions = random.sample(range(len(genome)), n_snps)

    for pos in positions:
        original = genome[pos]
        # 원래 염기를 제외한 3개 중 하나로 치환 → 반드시 변이가 되도록 보장한다.
        genome[pos] = random.choice(list(bases - {original}))

    return "".join(genome)


def save_sample(genome, out_path):
    """sample 게놈을 FASTA 형식(80자 줄바꿈)으로 저장."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write(">sample\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + "\n")
    print(f"Saved: {out_path} ({len(genome):,} bp)")


def generate_reads(genome, m, l, seed=None):
    """sample 게놈에서 길이 l짜리 read를 m개 무작위 위치에서 떼어낸다.

    Returns: (시작위치, 서열) 튜플의 리스트. 시작위치는 truth 정답으로 사용된다.
    """
    rng = random.Random(seed)
    reads = []
    max_start = len(genome) - l  # read가 게놈 끝을 넘지 않도록 시작 위치 상한

    for _ in range(m):
        pos = rng.randint(0, max_start)
        reads.append((pos, genome[pos:pos + l]))

    return reads


def save_reads(reads, out_path):
    """read 서열만 한 줄에 하나씩 저장 (위치 정보는 truth에 별도 저장)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        for _, seq in reads:
            f.write(seq + "\n")
    print(f"Saved: {out_path} ({len(reads):,} reads)")


def save_truth(reads, out_path):
    """각 read의 정답 시작 위치를 TSV로 저장 (정렬률·재현율 평가의 기준).

    read 이름은 save_reads의 줄 순서(idx)와 일치하도록 read_{idx:06d}로 부여한다.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("read_name\tref_start\n")
        for idx, (pos, _) in enumerate(reads):
            f.write(f"read_{idx:06d}\t{pos}\n")
    print(f"Truth:  {out_path} ({len(reads):,} entries)")


def load_genome(path):
    """FASTA 파일에서 헤더(>)를 제외한 서열만 이어붙여 로드."""
    sequence = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip())
    return "".join(sequence)


def save_reference(genome, path, header):
    """reference 게놈을 FASTA 형식(80자 줄바꿈)으로 저장."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + "\n")
    print(f"Saved: {path} ({len(genome):,} bp)")


def main():
    random.seed(SEED)

    # genome 크기 × SNP rate × coverage 모든 조합을 순회하며 데이터셋을 생성한다.
    # 이미 파일이 있으면 건너뛰어(skipping) 재실행 시 중복 생성을 막는다.
    for ds in DATASETS:
        n, label = ds["n"], ds["label"]
        ds_dir = os.path.join(DATA_DIR, label)
        ref_path = os.path.join(ds_dir, f"reference_{label}.festa")

        # 1. reference: ACGT를 무작위로 n개 생성
        if not os.path.exists(ref_path):
            print(f"\n[{label}] Generating random reference ({n:,} bp)...")
            genome = "".join(random.choices("ACGT", k=n))
            save_reference(genome, ref_path, header=f"reference_random_{label}")
        else:
            print(f"\n[{label}] Reference already exists, skipping.")

        reference = load_genome(ref_path)

        for snp in SNP_RATES:
            snp_dir = os.path.join(ds_dir, snp["label"])
            sample_path = os.path.join(snp_dir, "sample.festa")

            # 2. sample: reference에 SNP 주입 (정답 게놈)
            if not os.path.exists(sample_path):
                print(f"  [{snp['label']}] Generating sample (rate={snp['rate']})...")
                sample = inject_snps(reference, snp_rate=snp["rate"], seed=SEED)
                save_sample(sample, sample_path)
            else:
                print(f"  [{snp['label']}] Sample already exists, skipping.")
                sample = load_genome(sample_path)

            for coverage in COVERAGES:
                # read 개수 m = (genome 길이 / read 길이) × coverage
                m = n // L * coverage
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                truth_path = os.path.join(snp_dir, f"truth_{m}.tsv")
                # 3. reads + truth: sample에서 read를 떼어내고 정답 위치를 기록
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
