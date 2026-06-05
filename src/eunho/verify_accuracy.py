"""정확성 검증 — 이론 가능치 vs 실측 정렬률.

이론 가능치(Mappable) = 각 read를 truth(정답) 위치의 reference와 직접 비교해
mismatch ≤ D 인 read의 비율. = "이론상 매핑 가능한 read 비율".

Pigeonhole 알고리즘은 mismatch ≤ D 매핑을 빠짐없이 찾고(누락 0), random에선
오매핑이 없으므로(false 0), **정렬률(실측) == 이론 가능치(Mappable)** 가 성립해야 한다.
둘이 1:1로 일치하면 → 알고리즘이 이론 한계를 정확히 달성함을 입증.

운영점: SNP-matched D (1%→D=1, 3%→D=3, 5%→D=5).

실행(서버): python3 -m src.eunho.verify_accuracy
출력: results_random/mappable.tsv  (tag N SNP Coverage D Mappable AlignRate Diff)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.eunho.main import parse_fasta

GENOME_SIZES = [200_000, 600_000, 1_800_000]
COVERAGES    = [10, 20]
SNP_RATES    = [0.01, 0.03, 0.05]
SNP_TO_D     = {0.01: 1, 0.03: 3, 0.05: 5}   # SNP-matched 운영점
L            = 100
DATA_DIR     = "data/random"
RESULTS_DIR  = "results_random"

SIZE_LABEL = {200_000: "200K", 600_000: "600K", 1_800_000: "1800K"}
SNP_LABEL  = {0.01: "snp_1",   0.03: "snp_3",   0.05: "snp_5"}


def count_mm_le(a, b, D):
    """mismatch 수를 세되 D 초과 시 조기 종료."""
    mm = 0
    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > D:
                return mm
    return mm


def load_truth(path):
    truth = {}
    with open(path) as f:
        next(f)
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) >= 2:
                truth[p[0]] = int(p[1])
    return truth


def stream_reads(path):
    with open(path) as f:
        for idx, line in enumerate(f):
            seq = line.rstrip().upper()
            if seq:
                yield f"read_{idx:06d}", seq


def mappable_fraction(reference, reads_path, truth_path, D):
    """truth 위치에서 mismatch ≤ D 인 read 비율 (이론 가능치)."""
    truth = load_truth(truth_path)
    total = mappable = 0
    N = len(reference)
    for name, seq in stream_reads(reads_path):
        total += 1
        start = truth.get(name)
        if start is None or start < 0 or start + len(seq) > N:
            continue
        ref_sub = reference[start:start + len(seq)]
        if count_mm_le(ref_sub, seq, D) <= D:
            mappable += 1
    return mappable, total


def load_alignrate():
    """기존 summary.tsv 의 실측 정렬률(있으면) — 비교용."""
    path = os.path.join(RESULTS_DIR, "summary.tsv")
    rates = {}
    if not os.path.exists(path):
        return rates
    import csv
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("N"):
                rates[(r["N"], r["SNP"], r["Coverage"])] = r.get("AlignRate", "")
    return rates


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    align = load_alignrate()
    rows = []

    print(f"{'N':>6} {'SNP':>5} {'Cov':>5} {'D':>2} {'이론가능치':>10} {'실측정렬률':>10} {'차이':>7}")
    print("-" * 56)

    for Nv in GENOME_SIZES:
        size = SIZE_LABEL[Nv]
        ref_path = os.path.join(DATA_DIR, size, f"reference_{size}.festa")
        if not os.path.exists(ref_path):
            print(f"[SKIP] reference 없음: {ref_path}")
            continue
        reference = next(iter(parse_fasta(ref_path).values()))

        for snp in SNP_RATES:
            d = SNP_TO_D[snp]
            snp_dir = os.path.join(DATA_DIR, size, SNP_LABEL[snp])
            for cov in COVERAGES:
                m = Nv * cov // 100
                reads_path = os.path.join(snp_dir, f"reads_{m}.txt")
                truth_path = os.path.join(snp_dir, f"truth_{m}.tsv")
                if not (os.path.exists(reads_path) and os.path.exists(truth_path)):
                    print(f"[SKIP] reads/truth 없음: {reads_path}")
                    continue

                mappable, total = mappable_fraction(reference, reads_path, truth_path, d)
                mp = mappable / total * 100 if total else 0.0
                snp_lbl = f"{int(snp*100)}%"; cov_lbl = f"{cov}x"
                al = align.get((size, snp_lbl, cov_lbl), "")
                diff = ""
                if al.endswith("%"):
                    diff = f"{mp - float(al.rstrip('%')):+.2f}"
                rows.append((size, snp_lbl, cov_lbl, d, f"{mp:.2f}%", al, diff))
                print(f"{size:>6} {snp_lbl:>5} {cov_lbl:>5} {d:>2} "
                      f"{mp:>9.2f}% {al:>10} {diff:>7}")

    out = os.path.join(RESULTS_DIR, "mappable.tsv")
    with open(out, "w") as f:
        f.write("tag\tN\tSNP\tCoverage\tD\tMappable\tAlignRate\tDiff\n")
        for size, snp_lbl, cov_lbl, d, mp, al, diff in rows:
            f.write(f"random\t{size}\t{snp_lbl}\t{cov_lbl}\t{d}\t{mp}\t{al}\t{diff}\n")
    print(f"\n→ {out}")
    print("기대: 이론 가능치 == 실측 정렬률 (차이 ≈ 0) → 누락 0·오매핑 0 입증")


if __name__ == "__main__":
    main()
