import time
import os

from .main import run_mapping
from .evaluate import load_truth, load_result, evaluate

# ── 실험 설정 ────────────────────────────────────────────
GENOME_SIZES = [200_000, 600_000, 1_800_000]
COVERAGES    = [10, 20]
SNP_RATES    = [0.01, 0.03, 0.05]
DATA_DIR     = "data"
RESULTS_DIR  = "results"
# ─────────────────────────────────────────────────────────

SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL  = {0.01: "snp_1",   0.03: "snp_3",   0.05: "snp_5"}


def reads_count(N, cov, L=100):
    return N * cov // L


def main():
    results = []

    for N in GENOME_SIZES:
        size   = SIZE_LABEL[N]
        ref_path = os.path.join(DATA_DIR, size, f"reference_{size}.festa")

        for snp in SNP_RATES:
            snp_dir = os.path.join(DATA_DIR, size, SNP_LABEL[snp])

            for cov in COVERAGES:
                m          = reads_count(N, cov)
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                truth_path = os.path.join(snp_dir, f"truth_{m}.tsv")

                out_dir    = os.path.join(RESULTS_DIR, f"{size}_{SNP_LABEL[snp]}_cov{cov}x")
                os.makedirs(out_dir, exist_ok=True)
                result_path = os.path.join(out_dir, "result.tsv")

                label = f"N={N//1000}K  SNP={int(snp*100)}%  {cov}x"
                print(f"[실행] {label}", flush=True)

                if not os.path.exists(reads_path):
                    print(f"       [SKIP] reads 없음: {reads_path}")
                    continue

                # 정렬
                t0 = time.time()
                run_mapping(ref_path, reads_path, result_path)
                elapsed = time.time() - t0

                # result.tsv에서 정렬률 직접 계산
                mapped, total = 0, 0
                with open(result_path) as f:
                    next(f)  # 헤더 스킵
                    for line in f:
                        total += 1
                        if line.split("\t")[1] != "*":
                            mapped += 1

                entry = {
                    "N":      size,
                    "SNP":    f"{int(snp*100)}%",
                    "cov":    f"{cov}x",
                    "시간":   f"{elapsed:.1f}s",
                    "정렬률": f"{mapped/total*100:.1f}%" if total else "N/A",
                    "정밀도": "N/A",
                    "재현율": "N/A",
                }

                # 정밀도·재현율 (truth 파일이 있을 때만)
                if os.path.exists(truth_path):
                    metrics = evaluate(load_truth(truth_path), load_result(result_path))
                    entry["정밀도"] = f"{metrics['precision']*100:.1f}%"
                    entry["재현율"] = f"{metrics['recall']*100:.1f}%"

                results.append(entry)
                print(f"       정렬률={entry['정렬률']}  정밀도={entry['정밀도']}  재현율={entry['재현율']}  ({elapsed:.1f}s)")

    # 결과 테이블 출력
    print("\n" + "="*70)
    print(f"{'N':>8} {'SNP':>5} {'Coverage':>8} {'정렬률':>8} {'정밀도':>8} {'재현율':>8} {'시간':>8}")
    print("-"*70)
    for r in results:
        print(f"{r['N']:>8} {r['SNP']:>5} {r['cov']:>8} {r['정렬률']:>8} {r['정밀도']:>8} {r['재현율']:>8} {r['시간']:>8}")


if __name__ == "__main__":
    main()
