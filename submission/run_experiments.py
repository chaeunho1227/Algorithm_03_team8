import time
import os

from main import run_mapping, parse_fasta
from restoration import load_reads_by_name, restoration_and_snp_recall

DATA_DIR    = "data"
RESULTS_DIR = "results"

GENOME_SIZES = [200_000, 600_000, 1_800_000]
COVERAGES    = [10, 20]
SNP_RATES    = [0.01, 0.03, 0.05]
SNP_TO_D     = {0.01: 1, 0.03: 3, 0.05: 5}

SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL  = {0.01: "snp_1", 0.03: "snp_3", 0.05: "snp_5"}


def run_dataset(data_dir, results_dir):
    results = []
    for N in GENOME_SIZES:
        size = SIZE_LABEL[N]
        ref_path = os.path.join(data_dir, size, f"reference_{size}.festa")
        reference = next(iter(parse_fasta(ref_path).values()))

        for snp in SNP_RATES:
            d = SNP_TO_D[snp]
            snp_dir = os.path.join(data_dir, size, SNP_LABEL[snp])
            sample_path = os.path.join(snp_dir, "sample.festa")
            sample = next(iter(parse_fasta(sample_path).values())) if os.path.exists(sample_path) else None

            for cov in COVERAGES:
                m = N // 100 * cov
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                if not os.path.exists(reads_path):
                    print(f"[SKIP] {reads_path}")
                    continue

                out_dir = os.path.join(results_dir, f"{size}_{SNP_LABEL[snp]}_cov{cov}x")
                os.makedirs(out_dir, exist_ok=True)
                result_path = os.path.join(out_dir, "result.tsv")

                print(f"\n[실행] N={N//1000}K  SNP={int(snp*100)}%  D={d}  {cov}x")
                t0 = time.time()
                run_mapping(ref_path, reads_path, result_path, d=d)
                elapsed = time.time() - t0

                mapped = total = 0
                with open(result_path) as f:
                    next(f)
                    for line in f:
                        total += 1
                        if line.split("\t")[1] != "*":
                            mapped += 1

                entry = {
                    "N": size, "SNP": f"{int(snp*100)}%", "cov": f"{cov}x",
                    "정렬률": f"{mapped/total*100:.1f}%" if total else "N/A",
                    "재현율": "N/A", "복원정확도": "N/A",
                    "시간": f"{elapsed:.1f}s",
                }
                if sample is not None:
                    reads = load_reads_by_name(reads_path)
                    acc, rec = restoration_and_snp_recall(reference, sample, result_path, reads)
                    entry["재현율"] = f"{rec*100:.2f}%"
                    entry["복원정확도"] = f"{acc*100:.2f}%"

                results.append(entry)
                print(f"정렬률={entry['정렬률']}  재현율={entry['재현율']}  복원정확도={entry['복원정확도']}  시간={entry['시간']}")

    return results


def write_summary(results_dir, results):
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "summary.tsv"), "w") as f:
        f.write("N\tSNP\tCoverage\tAlignRate\tSNPRecall\tRestoreAcc\tTime\n")
        for r in results:
            f.write("\t".join([r["N"], r["SNP"], r["cov"], r["정렬률"], r["재현율"], r["복원정확도"], r["시간"]]) + "\n")
    print(f"요약 저장: {os.path.join(results_dir, 'summary.tsv')}")


def main():
    all_results = run_dataset(DATA_DIR, RESULTS_DIR)
    write_summary(RESULTS_DIR, all_results)
    print(f"\n{'N':>6} {'SNP':>5} {'Cov':>6} {'정렬률':>8} {'재현율':>8} {'복원정확도':>10} {'시간':>8}")
    print("-" * 60)
    for r in all_results:
        print(f"{r['N']:>6} {r['SNP']:>5} {r['cov']:>6} {r['정렬률']:>8} {r['재현율']:>8} {r['복원정확도']:>10} {r['시간']:>8}")


if __name__ == "__main__":
    main()
