"""SNP 단위 평가 — Voting 정책으로 매핑 결과에서 SNP를 콜한 뒤 truth와 비교.

전제:
- 기존 매핑 결과 TSV (submission/results/.../result.tsv)가 이미 존재.
- 재매핑 없이 SNP calling + evaluation만 수행.
"""
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from main import parse_fasta, stream_reads


# ======== 설정 ========
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
OUT_DIR = os.path.join(os.path.dirname(__file__), "results_snp")
os.makedirs(OUT_DIR, exist_ok=True)

GENOME_SIZES = [200_000, 600_000, 1_800_000]
COVERAGES = [10, 20]
SNP_RATES = [0.01, 0.03, 0.05]
SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL = {0.01: "snp_1", 0.03: "snp_3", 0.05: "snp_5"}


# ======== SNP truth 생성 ========
def generate_truth_snps(reference, sample):
    """Reference vs sample 직접 비교 — SNP 위치 집합 반환."""
    truth = set()
    n = min(len(reference), len(sample))
    for i in range(n):
        if reference[i] != sample[i]:
            truth.add(i)
    return truth


# ======== SNP calling: Voting 정책 ========
def call_snps(result_path, reads_path, reference):
    """매핑 결과로부터 voting 기반 SNP 콜.

    각 reference 위치 p에서 mismatch read 수가 depth/2를 초과하면 SNP로 콜.
    """
    # reads를 dict로 미리 로드 (read_name → seq)
    reads_dict = {}
    L = 0
    for name, seq in stream_reads(reads_path):
        reads_dict[name] = seq
        if L == 0:
            L = len(seq)

    depth = defaultdict(int)
    mm_count = defaultdict(int)

    # 매핑 결과 순회
    with open(result_path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            read_name, chrom, pos_str = parts[0], parts[1], parts[2]
            if chrom == "*":
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            if pos < 0 or pos + L > len(reference):
                continue
            read = reads_dict.get(read_name)
            if read is None:
                continue

            # 위치별 depth + mismatch 집계
            for j in range(L):
                p = pos + j
                depth[p] += 1
                if reference[p] != read[j]:
                    mm_count[p] += 1

    # Voting: mm_count[p] > depth[p] / 2
    called = set()
    for p, d in depth.items():
        if mm_count[p] * 2 > d:
            called.add(p)
    return called, depth, mm_count


# ======== Evaluation ========
def evaluate(called, truth):
    tp = len(called & truth)
    fp = len(called - truth)
    fn = len(truth - called)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return tp, fp, fn, precision, recall, f1


# ======== 메인 ========
def main():
    rows = []

    for N in GENOME_SIZES:
        size = SIZE_LABEL[N]
        ref_path = os.path.join(DATA_DIR, size, f"reference_{size}.festa")
        reference = next(iter(parse_fasta(ref_path).values()))

        for snp_rate in SNP_RATES:
            snp_lbl = SNP_LABEL[snp_rate]
            sample_path = os.path.join(DATA_DIR, size, snp_lbl, "sample.festa")
            sample = next(iter(parse_fasta(sample_path).values()))

            truth_snps = generate_truth_snps(reference, sample)

            for cov in COVERAGES:
                m = N * cov // 100
                reads_path = os.path.join(DATA_DIR, size, snp_lbl, f"reads_{m}.txt")
                result_path = os.path.join(RESULTS_DIR, f"{size}_{snp_lbl}_cov{cov}x", "result.tsv")

                if not os.path.exists(result_path):
                    print(f"[SKIP] {result_path} 없음")
                    continue

                label = f"N={size}  SNP={int(snp_rate*100)}%  cov={cov}x"
                print(f"\n=== {label} ===  truth_SNPs={len(truth_snps):,}")

                t0 = time.time()
                called, depth, _ = call_snps(result_path, reads_path, reference)
                t_call = time.time() - t0

                tp, fp, fn, prec, rec, f1 = evaluate(called, truth_snps)
                # SNP 위치 중 매핑이 한 번도 안 깔린 비율 (이론적 상한 계산용)
                uncovered_truth = sum(1 for s in truth_snps if depth.get(s, 0) == 0)

                row = {
                    "Genome": size,
                    "SNP": f"{int(snp_rate*100)}%",
                    "Coverage": f"{cov}x",
                    "truth": len(truth_snps),
                    "called": len(called),
                    "TP": tp,
                    "FP": fp,
                    "FN": fn,
                    "uncovered_truth": uncovered_truth,
                    "precision": f"{prec*100:.2f}%",
                    "recall": f"{rec*100:.2f}%",
                    "F1": f"{f1*100:.2f}%",
                    "calling_time": f"{t_call:.2f}s",
                }
                rows.append(row)
                print(f"  truth={len(truth_snps):,}  called={len(called):,}  TP={tp:,}  FP={fp}  FN={fn:,}")
                print(f"  precision={row['precision']}  recall={row['recall']}  F1={row['F1']}")
                print(f"  (uncovered truth={uncovered_truth},  calling time={t_call:.2f}s)")

    # ── 결과 저장 ──
    out_tsv = os.path.join(OUT_DIR, "snp_summary.tsv")
    with open(out_tsv, "w") as f:
        f.write("\t".join(rows[0].keys()) + "\n")
        for r in rows:
            f.write("\t".join(str(v) for v in r.values()) + "\n")
    print(f"\n✔ 저장: {out_tsv}")

    # ── 요약 표 ──
    print("\n" + "=" * 100)
    print(f"{'Genome':<8}{'SNP':<6}{'Cov':<6}{'Truth':<8}{'Called':<8}{'TP':<8}{'FP':<6}{'FN':<8}{'Prec':<10}{'Recall':<10}{'F1':<10}")
    print("-" * 100)
    for r in rows:
        print(f"{r['Genome']:<8}{r['SNP']:<6}{r['Coverage']:<6}{r['truth']:<8}{r['called']:<8}{r['TP']:<8}{r['FP']:<6}{r['FN']:<8}{r['precision']:<10}{r['recall']:<10}{r['F1']:<10}")


if __name__ == "__main__":
    main()
