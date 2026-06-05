"""정확도 회복 곡선 측정 (200K · 10x).

D=1/3/5/7/9 × SNP 1%/3%/5% = 15개 매핑을 새로 수행한 뒤,
각 조합의 정확도 = (전체 게놈 − uncovered SNP) / 전체 게놈 을 산출한다.

실행: python3 accuracy_recovery.py
출력:
- results/accuracy_recovery.tsv  (SNP × D 정확도 매트릭스)
- results/accuracy_recovery_full.tsv  (참고용 — 정확도/SNP recall 함께)
"""
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from main import parse_fasta, run_mapping


HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
OUT_DIR = os.path.join(HERE, "results")
TMP_DIR = os.path.join(HERE, "results", "_accuracy_recovery_tmp")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

# 측정 조건 — PPT 페이지 5와 동일
N_LABEL = "200K"
M = 20_000        # 200K · 10x
L = 100
D_VALUES = [1, 3, 5, 7, 9]
SNP_VALUES = [("snp_1", 0.01), ("snp_3", 0.03), ("snp_5", 0.05)]


def compute_depth(result_path, L_, N_):
    depth = defaultdict(int)
    with open(result_path) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[1]
            if chrom == "*":
                continue
            try:
                pos = int(parts[2])
            except ValueError:
                continue
            if pos < 0 or pos + L_ > N_:
                continue
            for j in range(L_):
                depth[pos + j] += 1
    return depth


def main():
    ref_path = os.path.join(DATA_DIR, N_LABEL, f"reference_{N_LABEL}.festa")
    reference = next(iter(parse_fasta(ref_path).values()))
    N = len(reference)
    print(f"Reference: {N:,} bp")

    results = {}  # snp_rate -> {D -> (accuracy, snp_recall)}

    for snp_lbl, snp_rate in SNP_VALUES:
        sample_path = os.path.join(DATA_DIR, N_LABEL, snp_lbl, "sample.festa")
        sample = next(iter(parse_fasta(sample_path).values()))
        truth_snps = {i for i in range(N) if reference[i] != sample[i]}
        n_truth = len(truth_snps)
        reads_path = os.path.join(DATA_DIR, N_LABEL, snp_lbl, f"reads_{M}.txt")
        print(f"\n=== SNP {int(snp_rate*100)}% (truth SNPs = {n_truth:,}) ===")

        results[snp_rate] = {}
        for D in D_VALUES:
            result_path = os.path.join(TMP_DIR, f"{snp_lbl}_D{D}.tsv")
            t0 = time.time()
            run_mapping(ref_path, reads_path, result_path, d=D)
            t_map = time.time() - t0

            depth = compute_depth(result_path, L, N)
            uncovered = sum(1 for s in truth_snps if depth.get(s, 0) == 0)
            accuracy = (N - uncovered) / N * 100
            snp_recall = (n_truth - uncovered) / n_truth * 100
            results[snp_rate][D] = (accuracy, snp_recall)
            print(f"  D={D}: 정확도={accuracy:.4f}%  SNP recall={snp_recall:.2f}%  (map {t_map:.2f}s)")

    # ── 매트릭스 형태 저장 (보고서용) ──
    out_matrix = os.path.join(OUT_DIR, "accuracy_recovery.tsv")
    with open(out_matrix, "w") as f:
        header = ["SNP/D"] + [f"D={D}" for D in D_VALUES]
        f.write("\t".join(header) + "\n")
        for _, snp_rate in SNP_VALUES:
            row = [f"{int(snp_rate*100)}%"]
            row += [f"{results[snp_rate][D][0]:.4f}" for D in D_VALUES]
            f.write("\t".join(row) + "\n")
    print(f"\n✔ 저장: {out_matrix}")

    # ── 전체 저장 (정확도 + SNP recall) ──
    out_full = os.path.join(OUT_DIR, "accuracy_recovery_full.tsv")
    with open(out_full, "w") as f:
        f.write("SNP\tD\t정확도(%)\tSNP_recall(%)\n")
        for _, snp_rate in SNP_VALUES:
            for D in D_VALUES:
                acc, rec = results[snp_rate][D]
                f.write(f"{int(snp_rate*100)}%\t{D}\t{acc:.4f}\t{rec:.2f}\n")
    print(f"✔ 저장: {out_full}")

    # ── 요약 표 출력 ──
    print("\n=== 정확도 회복 매트릭스 (200K · 10x) ===")
    header = ["SNP\\D"] + [f"D={D}" for D in D_VALUES]
    print("  " + "  ".join(f"{h:>9}" for h in header))
    for _, snp_rate in SNP_VALUES:
        row = [f"{int(snp_rate*100)}%"]
        row += [f"{results[snp_rate][D][0]:.4f}%" for D in D_VALUES]
        print("  " + "  ".join(f"{c:>9}" for c in row))


if __name__ == "__main__":
    main()
