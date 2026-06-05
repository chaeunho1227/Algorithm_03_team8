import time
import os

from main import run_mapping

#====== 실험 설정 ==============================================================

# reads / reference 파일이 들어 있는 루트 디렉토리
DATA_DIR = "data"

# 매핑 결과(result.tsv)를 저장할 루트 디렉토리
RESULTS_DIR = "results"

# 테스트할 genome 크기 (200_000 / 600_000 / 1_800_000)
GENOME_SIZES = [200_000, 600_000, 1_800_000]

# 테스트할 coverage (10 / 20)
COVERAGES = [10, 20]

# 테스트할 SNP rate (0.01 / 0.03 / 0.05)
SNP_RATES = [0.01, 0.03, 0.05]

# 반복 서열 필터 임계값 (0 = 필터 없음)
MAX_REPEAT = 500

# SNP rate별 허용 mismatch 수
# 고정: {0.01: 3, 0.03: 3, 0.05: 3}
# SNP-matched: {0.01: 1, 0.03: 3, 0.05: 5}
SNP_TO_D = {0.01: 1, 0.03: 3, 0.05: 5}

#====================================================================

SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL = {0.01: "snp_1", 0.03: "snp_3", 0.05: "snp_5"}


def load_truth(path):
    truth = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            truth[parts[0]] = int(parts[1])
    return truth


def load_result(path):
    result = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            result[parts[0]] = (parts[1], int(parts[2]))
    return result


def evaluate(truth, result):
    total = len(truth)
    mapped = 0
    correct = 0

    for name, true_start in truth.items():
        chrom, pred_start = result[name]
        if chrom == "*":
            continue
        mapped += 1
        if pred_start == true_start:
            correct += 1

    precision = correct / mapped
    recall = correct / total

    return {"precision": precision, "recall": recall}


def run_dataset(data_dir, results_dir):
    results = []

    for N in GENOME_SIZES:
        size = SIZE_LABEL[N]
        ref_path = os.path.join(data_dir, size, f"reference_{size}.festa")

        for snp in SNP_RATES:
            d = SNP_TO_D[snp]
            snp_dir = os.path.join(data_dir, size, SNP_LABEL[snp])

            for cov in COVERAGES:
                m = N * cov // 100
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                truth_path = os.path.join(snp_dir, f"truth_{m}.tsv")

                out_dir = os.path.join(results_dir, f"{size}_{SNP_LABEL[snp]}_cov{cov}x")
                os.makedirs(out_dir, exist_ok=True)
                result_path = os.path.join(out_dir, "result.tsv")

                label = f"N={N//1000}K  SNP={int(snp*100)}%  D={d}  {cov}x"
                print(f"\n===============================================================\n", flush=True)
                print(f"[실행] {label}")

                t0 = time.time()
                run_mapping(ref_path, reads_path, result_path, max_repeat=MAX_REPEAT, d=d)
                elapsed = time.time() - t0

                mapped, total = 0, 0
                with open(result_path) as f:
                    next(f)
                    for line in f:
                        total += 1
                        if line.split("\t")[1] != "*":
                            mapped += 1

                metrics = evaluate(load_truth(truth_path), load_result(result_path))
                entry = {
                    "N": size,
                    "SNP": f"{int(snp*100)}%",
                    "cov": f"{cov}x",
                    "시간": f"{elapsed:.1f}s",
                    "정렬률": f"{mapped/total*100:.1f}%",
                    "정밀도": f"{metrics['precision']*100:.1f}%",
                    "재현율": f"{metrics['recall']*100:.1f}%",
                }

                results.append(entry)
                print(f"정렬률={entry['정렬률']}  정밀도={entry['정밀도']}  재현율={entry['재현율']}  ({elapsed:.1f}s)")

    return results


def main():
    print(f"\n데이터셋: data={DATA_DIR}, results={RESULTS_DIR}")
    all_results = run_dataset(DATA_DIR, RESULTS_DIR)

    print()
    print(f"{'N':>6} {'SNP':>5} {'Coverage':>8} {'정렬률':>8} {'정밀도':>8} {'재현율':>8} {'시간':>8}")
    for r in all_results:
        print(f"{r['N']:>6} {r['SNP']:>5} {r['cov']:>8} "
              f"{r['정렬률']:>8} {r['정밀도']:>8} {r['재현율']:>8} {r['시간']:>8}")


if __name__ == "__main__":
    main()
