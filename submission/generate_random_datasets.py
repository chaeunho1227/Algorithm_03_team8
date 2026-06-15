import os
import random

DATA_DIR = "data"
L = 100
SEED = 1212

GENOME_SIZES = [(200_000, "200K"), (600_000, "600K"), (1_800_000, "1800K")]
SNP_RATES    = [(0.01, "snp_1"), (0.03, "snp_3"), (0.05, "snp_5")]
COVERAGES    = [10, 20]


def save_fasta(genome, path, header):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + "\n")
    print(f"Saved: {path} ({len(genome):,} bp)")


def load_genome(path):
    with open(path) as f:
        return "".join(l.strip() for l in f if not l.startswith(">"))


def inject_snps(genome, snp_rate):
    random.seed(SEED)
    genome = list(genome)
    bases = {"A", "C", "G", "T"}
    for pos in random.sample(range(len(genome)), int(len(genome) * snp_rate)):
        genome[pos] = random.choice(list(bases - {genome[pos]}))
    return "".join(genome)


def generate_reads(genome, m):
    rng = random.Random(SEED)
    max_start = len(genome) - L
    reads = []
    for _ in range(m):
        pos = rng.randint(0, max_start)
        reads.append((pos, genome[pos:pos + L]))
    return reads


def save_reads(reads, reads_path, truth_path):
    os.makedirs(os.path.dirname(reads_path), exist_ok=True)
    with open(reads_path, "w") as rf, open(truth_path, "w") as tf:
        tf.write("read_name\tref_start\n")
        for idx, (pos, seq) in enumerate(reads):
            rf.write(seq + "\n")
            tf.write(f"read_{idx:06d}\t{pos}\n")
    print(f"Saved: {reads_path} ({len(reads):,} reads)")


def main():
    random.seed(SEED)
    for n, size in GENOME_SIZES:
        ref_path = os.path.join(DATA_DIR, size, f"reference_{size}.festa")
        if not os.path.exists(ref_path):
            save_fasta("".join(random.choices("ACGT", k=n)), ref_path, f"reference_random_{size}")
        reference = load_genome(ref_path)

        for rate, snp_label in SNP_RATES:
            snp_dir = os.path.join(DATA_DIR, size, snp_label)
            sample_path = os.path.join(snp_dir, "sample.festa")
            if not os.path.exists(sample_path):
                save_fasta(inject_snps(reference, rate), sample_path, "sample")
            sample = load_genome(sample_path)

            for cov in COVERAGES:
                m = n // L * cov
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                truth_path = os.path.join(snp_dir, f"truth_{m}.tsv")
                if not os.path.exists(reads_path):
                    save_reads(generate_reads(sample, m), reads_path, truth_path)

    print("All datasets ready.")


if __name__ == "__main__":
    main()
