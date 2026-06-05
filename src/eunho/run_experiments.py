import time
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.eunho.main import run_mapping, parse_fasta
from src.eunho.restoration import load_reads_by_name, restoration_and_snp_recall

# ==============================================================
#  실험 설정
# ==============================================================

# 테스트할 genome 크기 (200_000 / 600_000 / 1_800_000)
GENOME_SIZES = [200_000, 600_000, 1_800_000]

# 테스트할 coverage (10 / 20)
COVERAGES = [10, 20]

# 테스트할 SNP rate (0.01 / 0.03 / 0.05)
SNP_RATES = [0.01, 0.03, 0.05]

# 반복 서열 필터 임계값 (0 = 필터 없음)
MAX_REPEAT = 500

# SNP rate별 허용 mismatch 수 D (운영점)
#   SNP-matched D: SNP rate에 맞춰 동적 설정 (1%→1, 3%→3, 5%→5)  ← 최종 운영점
#   고정 D=3      : {0.01: 3, 0.03: 3, 0.05: 3}
SNP_TO_D = {0.01: 1, 0.03: 3, 0.05: 5}

# 테스트할 데이터셋
DATASETS = [
    # {"data_dir": "data",        "results_dir": "results",        "tag": "chr1"},
    {"data_dir": "data/random", "results_dir": "results_random", "tag": "random"},
]

# ==============================================================
#  측정 지표 (4지표)  ※ 정밀도는 random에서 항상 100%라 제외
#   1. 정렬률    = 정렬된 read 수 / 전체 read 수            (완전성)
#   2. 재현율    = 탐지된 SNP 위치 수 / 실제 SNP 위치 수    (SNP 단위 탐지)
#   3. 복원정확도 = 일치 염기 수 / 전체 염기 수             (최종 복원 품질)
#   4. 시간      = 매핑 실행 시간 (초)                      (성능)
# ==============================================================

SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL  = {0.01: "snp_1",   0.03: "snp_3",   0.05: "snp_5"}


def run_dataset(data_dir, results_dir, tag):
    results = []

    for N in GENOME_SIZES:
        size      = SIZE_LABEL[N]
        ref_path  = os.path.join(data_dir, size, f"reference_{size}.festa")
        reference = next(iter(parse_fasta(ref_path).values()))

        for snp in SNP_RATES:
            d           = SNP_TO_D[snp]
            snp_dir     = os.path.join(data_dir, size, SNP_LABEL[snp])
            sample_path = os.path.join(snp_dir, "sample.festa")
            sample = None
            if os.path.exists(sample_path):
                sample = next(iter(parse_fasta(sample_path).values()))

            for cov in COVERAGES:
                m          = N * cov // 100
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")

                out_dir     = os.path.join(results_dir, f"{size}_{SNP_LABEL[snp]}_cov{cov}x")
                os.makedirs(out_dir, exist_ok=True)
                result_path = os.path.join(out_dir, "result.tsv")
                time_path   = os.path.join(out_dir, "time.txt")

                label = f"[{tag}] N={N//1000}K  SNP={int(snp*100)}%  D={d}  {cov}x"
                print(f"[실행] {label}", flush=True)

                if not os.path.exists(reads_path):
                    print(f"       [SKIP] reads 없음: {reads_path}")
                    continue

                # ── 매핑 (시간 측정) ──
                t0 = time.time()
                run_mapping(ref_path, reads_path, result_path, max_repeat=MAX_REPEAT, d=d)
                elapsed = time.time() - t0
                with open(time_path, "w") as f:
                    f.write(f"{elapsed:.3f}\n")

                # ── 정렬률 ──
                mapped, total = 0, 0
                with open(result_path) as f:
                    next(f)
                    for line in f:
                        total += 1
                        if line.split("\t")[1] != "*":
                            mapped += 1

                entry = {
                    "tag":        tag,
                    "N":          size,
                    "SNP":        f"{int(snp*100)}%",
                    "cov":        f"{cov}x",
                    "정렬률":     f"{mapped/total*100:.1f}%" if total else "N/A",
                    "재현율":     "N/A",
                    "복원정확도": "N/A",
                    "시간":       f"{elapsed:.1f}s",
                }

                # ── 재현율(SNP) + 복원정확도 (consensus 복원) ──
                if sample is not None:
                    reads = load_reads_by_name(reads_path)
                    acc, rec = restoration_and_snp_recall(reference, sample, result_path, reads)
                    entry["재현율"]     = f"{rec*100:.2f}%"
                    entry["복원정확도"] = f"{acc*100:.2f}%"

                results.append(entry)
                print(f"       정렬률={entry['정렬률']}  재현율={entry['재현율']}  "
                      f"복원정확도={entry['복원정확도']}  시간={entry['시간']}")

    return results


def write_summary(results_dir, results):
    summary_path = os.path.join(results_dir, "summary.tsv")
    with open(summary_path, "w") as f:
        f.write("tag\tN\tSNP\tCoverage\tAlignRate\tSNPRecall\tRestoreAcc\tTime\n")
        for r in results:
            f.write("\t".join([
                r["tag"], r["N"], r["SNP"], r["cov"],
                r["정렬률"], r["재현율"], r["복원정확도"], r["시간"],
            ]) + "\n")
    print(f"  → {summary_path}")


def main():
    all_results = []

    for ds in DATASETS:
        print(f"\n{'='*70}")
        print(f"  데이터셋: {ds['tag']}  ({ds['data_dir']})  —  SNP-matched D {SNP_TO_D}")
        print(f"{'='*70}")
        ds_results = run_dataset(ds["data_dir"], ds["results_dir"], ds["tag"])
        write_summary(ds["results_dir"], ds_results)
        all_results += ds_results

    print("\n" + "=" * 86)
    print(f"{'태그':>8} {'N':>6} {'SNP':>5} {'Cov':>5} "
          f"{'정렬률':>8} {'재현율':>8} {'복원정확도':>10} {'시간':>8}")
    print("-" * 86)
    for r in all_results:
        print(f"{r['tag']:>8} {r['N']:>6} {r['SNP']:>5} {r['cov']:>5} "
              f"{r['정렬률']:>8} {r['재현율']:>8} {r['복원정확도']:>10} {r['시간']:>8}")


if __name__ == "__main__":
    main()
