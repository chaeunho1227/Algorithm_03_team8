"""5지표 평가 — 18개 조합에 대해 정확도·정렬률·정밀도(SNP)·재현율(SNP) 측정.

전제: run_experiments.py가 먼저 실행되어 results/{Genome}_{SNP}_cov{Cov}x/result.tsv가 존재해야 함.

지표 정의 (PPT 그림 #8 기준):
- 정확도   = (전체 게놈 길이 − uncovered SNP 수) / 전체 게놈 길이      (염기 단위)
- 정렬률   = 매핑 성공 read / 전체 read                                (read 단위)
- 정밀도   = TP / (TP + FP)                                            (SNP 단위)
- 재현율   = TP / (TP + FN)                                            (SNP 단위)

SNP calling 정책: 매핑된 read에서 reference와 다른 base를 보고한 위치 중
voting (mismatch read 수 > depth/2) 기준으로 SNP 콜.

실행: python3 snp_evaluation.py
출력: results_snp/snp_summary.tsv
"""
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from main import parse_fasta, stream_reads


HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
RESULTS_DIR = os.path.join(HERE, "results")
OUT_DIR = os.path.join(HERE, "results_snp")
os.makedirs(OUT_DIR, exist_ok=True)

GENOME_SIZES = [200_000, 600_000, 1_800_000]
COVERAGES = [10, 20]
SNP_RATES = [0.01, 0.03, 0.05]
SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL = {0.01: "snp_1", 0.03: "snp_3", 0.05: "snp_5"}


def load_truth_read(path):
    truth = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                truth[parts[0]] = int(parts[1])
    return truth


def compute_truth_snps(reference, sample):
    n = min(len(reference), len(sample))
    return {i for i in range(n) if reference[i] != sample[i]}


def evaluate_one(reference, sample, result_path, truth_read_path, reads_path):
    """한 조합의 5지표 중 정확도·정렬률·정밀도(SNP)·재현율(SNP) 계산."""
    N = len(reference)
    truth_snps = compute_truth_snps(reference, sample)
    n_truth = len(truth_snps)
    truth_read = load_truth_read(truth_read_path)

    # reads 로드 (read_name → seq)
    reads_dict = {}
    L = 0
    for name, seq in stream_reads(reads_path):
        reads_dict[name] = seq
        if L == 0:
            L = len(seq)

    total = 0
    mapped = 0
    correct = 0
    depth = defaultdict(int)
    mm_count = defaultdict(int)

    with open(result_path) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            total += 1
            read_name, chrom, pos_str = parts[0], parts[1], parts[2]
            if chrom == "*":
                continue
            mapped += 1
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            if pos == truth_read.get(read_name):
                correct += 1
            if pos < 0 or pos + L > N:
                continue
            read = reads_dict.get(read_name)
            if read is None:
                continue
            for j in range(L):
                p = pos + j
                depth[p] += 1
                if reference[p] != read[j]:
                    mm_count[p] += 1

    # SNP calling: voting (mismatch > depth/2)
    called_snps = {p for p, d in depth.items() if mm_count[p] * 2 > d}

    tp = len(called_snps & truth_snps)
    fp = len(called_snps - truth_snps)
    fn = len(truth_snps - called_snps)
    uncovered_snp = sum(1 for s in truth_snps if depth.get(s, 0) == 0)

    return {
        "total": total,
        "mapped": mapped,
        "correct": correct,
        "alignment_rate": mapped / total if total > 0 else 0.0,
        "read_precision": correct / mapped if mapped > 0 else 0.0,
        "snp_precision": tp / (tp + fp) if (tp + fp) > 0 else 0.0,
        "snp_recall": tp / n_truth if n_truth > 0 else 0.0,
        "accuracy": (N - uncovered_snp) / N if N > 0 else 0.0,
        "tp": tp, "fp": fp, "fn": fn,
        "n_truth": n_truth,
        "uncovered_snp": uncovered_snp,
    }


def main():
    rows = []
    print(f"평가 대상: {len(GENOME_SIZES) * len(SNP_RATES) * len(COVERAGES)}개 조합")

    for N in GENOME_SIZES:
        size = SIZE_LABEL[N]
        ref_path = os.path.join(DATA_DIR, size, f"reference_{size}.festa")
        reference = next(iter(parse_fasta(ref_path).values()))

        for snp_rate in SNP_RATES:
            snp_lbl = SNP_LABEL[snp_rate]
            sample_path = os.path.join(DATA_DIR, size, snp_lbl, "sample.festa")
            sample = next(iter(parse_fasta(sample_path).values()))

            for cov in COVERAGES:
                m = N * cov // 100
                reads_path = os.path.join(DATA_DIR, size, snp_lbl, f"reads_{m}.txt")
                truth_path = os.path.join(DATA_DIR, size, snp_lbl, f"truth_{m}.tsv")
                result_path = os.path.join(RESULTS_DIR, f"{size}_{snp_lbl}_cov{cov}x", "result.tsv")

                if not os.path.exists(result_path):
                    print(f"[SKIP] {result_path} 없음 — run_experiments.py를 먼저 실행하세요.")
                    continue

                t0 = time.time()
                m_ = evaluate_one(reference, sample, result_path, truth_path, reads_path)
                t_eval = time.time() - t0

                row = {
                    "Genome": size,
                    "SNP": f"{int(snp_rate*100)}%",
                    "Coverage": f"{cov}x",
                    "정확도": f"{m_['accuracy']*100:.4f}%",
                    "정렬률": f"{m_['alignment_rate']*100:.2f}%",
                    "정밀도(SNP)": f"{m_['snp_precision']*100:.2f}%",
                    "재현율(SNP)": f"{m_['snp_recall']*100:.2f}%",
                    "TP": m_["tp"],
                    "FP": m_["fp"],
                    "FN": m_["fn"],
                    "truth_SNPs": m_["n_truth"],
                    "uncovered_SNPs": m_["uncovered_snp"],
                    "평가시간": f"{t_eval:.2f}s",
                }
                rows.append(row)
                print(f"=== {size}·{int(snp_rate*100)}%·{cov}x ===  "
                      f"정확도={row['정확도']}  정렬률={row['정렬률']}  "
                      f"정밀도={row['정밀도(SNP)']}  재현율={row['재현율(SNP)']}  "
                      f"(eval {t_eval:.1f}s)")

    if not rows:
        print("평가할 결과가 없습니다. results/ 디렉토리에 result.tsv가 있는지 확인하세요.")
        return

    out_path = os.path.join(OUT_DIR, "snp_summary.tsv")
    with open(out_path, "w") as f:
        f.write("\t".join(rows[0].keys()) + "\n")
        for r in rows:
            f.write("\t".join(str(v) for v in r.values()) + "\n")
    print(f"\n✔ 저장: {out_path}")

    # 요약 표 출력
    print("\n" + "=" * 110)
    print(f"{'Genome':<8}{'SNP':<6}{'Cov':<6}{'정확도':<12}{'정렬률':<10}{'정밀도(SNP)':<14}{'재현율(SNP)':<14}")
    print("-" * 110)
    for r in rows:
        print(f"{r['Genome']:<8}{r['SNP']:<6}{r['Coverage']:<6}"
              f"{r['정확도']:<12}{r['정렬률']:<10}{r['정밀도(SNP)']:<14}{r['재현율(SNP)']:<14}")


if __name__ == "__main__":
    main()
