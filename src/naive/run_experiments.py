import time
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.naive.main import run_mapping
from src.eunho.evaluate import load_truth, load_result, evaluate

# ==============================================================
#  실험 설정
# ==============================================================

# 테스트할 genome 크기 (200_000 / 600_000 / 1_800_000)
GENOME_SIZES = [200_000, 600_000, 1_800_000]

# 테스트할 coverage (10 / 20)
COVERAGES = [10, 20]

# 테스트할 SNP rate (0.01 / 0.03 / 0.05)
SNP_RATES = [0.01, 0.03, 0.05]

# 테스트할 데이터셋
# data_dir:    reads/reference 파일 위치
# results_dir: 결과 저장 위치
# tag:         결과 테이블에 표시될 이름
DATASETS = [
    {"data_dir": "data",        "results_dir": "results_naive",        "tag": "chr1"},
    {"data_dir": "data/random", "results_dir": "results_naive_random", "tag": "random"},
]

# ==============================================================

SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL  = {0.01: "snp_1",   0.03: "snp_3",   0.05: "snp_5"}


def reads_count(N, cov, L=100):
    return N * cov // L


def run_dataset(data_dir, results_dir, tag):
    results = []

    for N in GENOME_SIZES:
        size     = SIZE_LABEL[N]
        ref_path = os.path.join(data_dir, size, f"reference_{size}.festa")

        for snp in SNP_RATES:
            snp_dir    = os.path.join(data_dir, size, SNP_LABEL[snp])

            for cov in COVERAGES:
                m          = reads_count(N, cov)
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                truth_path = os.path.join(snp_dir, f"truth_{m}.tsv")

                out_dir     = os.path.join(results_dir, f"{size}_{SNP_LABEL[snp]}_cov{cov}x")
                os.makedirs(out_dir, exist_ok=True)
                result_path = os.path.join(out_dir, "result.tsv")

                label = f"[{tag}] N={N//1000}K  SNP={int(snp*100)}%  {cov}x"
                print(f"[실행] {label}", flush=True)

                if not os.path.exists(reads_path):
                    print(f"       [SKIP] reads 없음: {reads_path}")
                    continue

                t0 = time.time()
                run_mapping(ref_path, reads_path, result_path)
                elapsed = time.time() - t0

                mapped, total = 0, 0
                with open(result_path) as f:
                    next(f)
                    for line in f:
                        total += 1
                        if line.split("\t")[1] != "*":
                            mapped += 1

                entry = {
                    "tag":    tag,
                    "N":      size,
                    "SNP":    f"{int(snp*100)}%",
                    "cov":    f"{cov}x",
                    "시간":   f"{elapsed:.1f}s",
                    "정렬률": f"{mapped/total*100:.1f}%" if total else "N/A",
                    "정밀도": "N/A",
                    "재현율": "N/A",
                }

                if os.path.exists(truth_path):
                    metrics = evaluate(load_truth(truth_path), load_result(result_path))
                    entry["정밀도"] = f"{metrics['precision']*100:.1f}%"
                    entry["재현율"] = f"{metrics['recall']*100:.1f}%"

                results.append(entry)
                print(f"       정렬률={entry['정렬률']}  정밀도={entry['정밀도']}  재현율={entry['재현율']}  ({elapsed:.1f}s)")

    return results


def main():
    for ds in DATASETS:
        print(f"\n{'='*70}")
        print(f"  데이터셋: {ds['tag']}  ({ds['data_dir']})")
        print(f"{'='*70}")
        entries = run_dataset(ds["data_dir"], ds["results_dir"], ds["tag"])

        print("\n" + "="*78)
        print(f"{'태그':>8} {'N':>6} {'SNP':>5} {'Coverage':>8} {'정렬률':>8} {'정밀도':>8} {'재현율':>8} {'시간':>8}")
        print("-"*78)
        for r in entries:
            print(f"{r['tag']:>8} {r['N']:>6} {r['SNP']:>5} {r['cov']:>8} "
                  f"{r['정렬률']:>8} {r['정밀도']:>8} {r['재현율']:>8} {r['시간']:>8}")


if __name__ == "__main__":
    main()
