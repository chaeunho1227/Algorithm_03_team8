import random
import os

# ── 실험 설정 ────────────────────────────────────────────
REF_PATH     = "data/ch1/chr1_chunk.festa"             # 원본 레퍼런스 (1.8M bp)
GENOME_SIZES = [200_000, 600_000, 1_800_000]       # N (앞에서 N bp 사용)
COVERAGES    = [10, 20]                            # x coverage
SNP_RATES    = [0.01, 0.03, 0.05]                 # 각 위치 독립 변이 확률
READ_LENGTH  = 100                                 # L (bp)
SEED         = 1227
OUT_DIR      = "data/experiments"
# ─────────────────────────────────────────────────────────


def load_genome(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def generate_reads(genome, num_reads, read_length, snp_rate, seed=SEED):
    rng = random.Random(seed)
    n = len(genome)
    reads = []
    truth = []

    for i in range(num_reads):
        name = f"read_{i:06d}"

        start    = rng.randint(0, n - read_length)
        original = genome[start:start + read_length]

        strand = rng.choice(["+", "-"])
        seq = original if strand == "+" else reverse_complement(original)

        seq = list(seq)
        for pos in range(read_length):
            if rng.random() < snp_rate:
                seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
        ref_seq = original if strand == "+" else reverse_complement(original)
        num_mut = sum(a != b for a, b in zip(seq, ref_seq))
        seq = "".join(seq)

        reads.append((name, seq))
        truth.append((name, start, strand, num_mut))

    return reads, truth


def write_fasta(reads, path):
    with open(path, "w") as f:
        for name, seq in reads:
            f.write(f">{name}\n{seq}\n")


def write_truth(truth, path):
    with open(path, "w") as f:
        f.write("read_name\tref_start\tstrand\tnum_mutations\n")
        for name, start, strand, num_mut in truth:
            f.write(f"{name}\t{start}\t{strand}\t{num_mut}\n")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"레퍼런스 로드: {REF_PATH}")
    full_genome = load_genome(REF_PATH)
    print(f"  전체 길이: {len(full_genome):,} bp\n")

    for N in GENOME_SIZES:
        genome  = full_genome[:N]
        n_label = f"N{N // 1000}K"

        for snp in SNP_RATES:
            snp_label = f"snp{int(snp * 100)}"

            max_reads = N * max(COVERAGES) // READ_LENGTH
            print(f"생성: {n_label}/{snp_label}  (N={N:,}, SNP={int(snp*100)}%, reads={max_reads:,})")

            reads, truth = generate_reads(
                genome,
                num_reads=max_reads,
                read_length=READ_LENGTH,
                snp_rate=snp,
                seed=SEED + 1,
            )

            for cov in COVERAGES:
                num_reads = N * cov // READ_LENGTH
                cov_label = f"cov{cov}x"

                # data/experiments/N200K/snp1/cov10x/
                out_dir = os.path.join(OUT_DIR, n_label, snp_label, cov_label)
                os.makedirs(out_dir, exist_ok=True)

                # 레퍼런스는 N별 공용
                ref_path = os.path.join(OUT_DIR, n_label, "ref.fasta")
                if not os.path.exists(ref_path):
                    os.makedirs(os.path.join(OUT_DIR, n_label), exist_ok=True)
                    with open(ref_path, "w") as f:
                        f.write(f">chr1\n{genome}\n")

                write_fasta(reads[:num_reads], os.path.join(out_dir, "reads.fasta"))
                write_truth(truth[:num_reads], os.path.join(out_dir, "truth.tsv"))

    print(f"\n완료 → {OUT_DIR}/")


if __name__ == "__main__":
    main()
