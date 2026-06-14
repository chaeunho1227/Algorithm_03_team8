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
    """매핑 결과로부터 reference 위치별 커버리지 깊이(depth)를 계산.

    각 매핑된 read가 덮는 [pos, pos+L) 구간의 모든 위치 카운트를 1씩 올린다.
    depth[i] == 0 이면 그 위치를 덮는 read가 하나도 없다는 뜻 → 복원 불가.

    Returns: position -> depth 딕셔너리 (커버된 위치만 보관).
    """
    depth = defaultdict(int)
    with open(result_path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[1]
            if chrom == "*":  # 매핑 실패한 read는 제외
                continue
            try:
                pos = int(parts[2])
            except ValueError:
                continue
            if pos < 0 or pos + L_ > N_:  # 게놈 범위를 벗어나는 매핑은 제외
                continue
            for j in range(L_):
                depth[pos + j] += 1
    return depth


def main():
    """SNP × D 조합별 정확도를 매핑→측정→저장하는 전체 파이프라인."""
    ref_path = os.path.join(DATA_DIR, N_LABEL, f"reference_{N_LABEL}.festa")
    reference = next(iter(parse_fasta(ref_path).values()))
    N = len(reference)
    print(f"Reference: {N:,} bp")

    results = {}  # snp_rate -> {D -> (accuracy, snp_recall)}

    for snp_lbl, snp_rate in SNP_VALUES:
        sample_path = os.path.join(DATA_DIR, N_LABEL, snp_lbl, "sample.festa")
        sample = next(iter(parse_fasta(sample_path).values()))
        # 실제 SNP 위치 = reference와 sample 염기가 다른 위치
        truth_snps = {i for i in range(N) if reference[i] != sample[i]}
        n_truth = len(truth_snps)
        reads_path = os.path.join(DATA_DIR, N_LABEL, snp_lbl, f"reads_{M}.txt")
        print(f"\n=== SNP {int(snp_rate*100)}% (truth SNPs = {n_truth:,}) ===")

        # D(허용 mismatch)를 1→9로 키워가며 정확도가 회복되는 곡선을 측정한다.
        results[snp_rate] = {}
        for D in D_VALUES:
            result_path = os.path.join(TMP_DIR, f"{snp_lbl}_D{D}.tsv")
            t0 = time.time()
            run_mapping(ref_path, reads_path, result_path, d=D)
            t_map = time.time() - t0

            depth = compute_depth(result_path, L, N)
            # SNP 위치를 덮는 read가 하나도 없으면(depth 0) 그 변이는 복원할 수 없다.
            uncovered = sum(1 for s in truth_snps if depth.get(s, 0) == 0)
            # 정확도   = 복원 가능한(=커버된) 전체 게놈 비율
            accuracy = (N - uncovered) / N * 100
            # SNP recall = 커버되어 탐지 가능한 SNP 비율
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
