import os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from src.eunho.main import parse_fasta
from src.eunho.restoration import load_reads_by_name, restoration_accuracy

SIZE_BP = {"200K": 200_000, "600K": 600_000, "1800K": 1_800_000}
PATTERN = re.compile(r"^(200K|600K|1800K)_(snp_[135])(?:_d\d+)?_cov(\d+)x$")

# 4지표 측정
#   정렬률    = 정렬된 read 수 / 전체 read 수
#   정밀도    = 제자리에 붙인 read 수 / 매핑한 read 수
#   복원정확도 = 일치 염기 수 / 전체 염기 수 (consensus 복원 후 sample 과 비교)
#   시간      = 매핑 실행 시간 (초)  ← results/{combo}/time.txt 가 있으면 표시


def analyze(results_dir):
    data_dir = "data/random" if "random" in results_dir else "data"
    tag = "random" if "random" in results_dir else "chr1"

    ref_cache = {}      # size -> reference
    sample_cache = {}   # (size, snp) -> sample

    entries = []
    for subdir in os.listdir(results_dir):
        m = PATTERN.match(subdir)
        result_path = os.path.join(results_dir, subdir, "result.tsv")
        if not m or not os.path.exists(result_path):
            continue

        size, snp, cov = m.group(1), m.group(2), int(m.group(3))
        reads = SIZE_BP[size] * cov // 100
        truth_path = os.path.join(data_dir, size, snp, f"truth_{reads}.tsv")
        reads_path = os.path.join(data_dir, size, snp, f"reads_{reads}.txt")
        ref_path = os.path.join(data_dir, size, f"reference_{size}.festa")
        sample_path = os.path.join(data_dir, size, snp, "sample.festa")

        # ── 정렬률 ──
        total = mapped = 0
        with open(result_path) as f:
            next(f)
            for line in f:
                total += 1
                if line.split("\t")[1] != "*":
                    mapped += 1

        # ── 정밀도 ──
        prec = "N/A"
        if os.path.exists(truth_path):
            truth = {}
            with open(truth_path) as f:
                next(f)
                for line in f:
                    p = line.rstrip().split("\t")
                    truth[p[0]] = int(p[1])
            correct = 0
            with open(result_path) as f:
                next(f)
                for line in f:
                    p = line.rstrip().split("\t")
                    if p[1] != "*" and truth.get(p[0]) == int(p[2]):
                        correct += 1
            prec = f"{correct/mapped*100:.1f}%" if mapped else "N/A"

        # ── 복원정확도 (consensus) ──
        restore = "N/A"
        if os.path.exists(ref_path) and os.path.exists(sample_path) and os.path.exists(reads_path):
            if size not in ref_cache:
                ref_cache[size] = next(iter(parse_fasta(ref_path).values()))
            key = (size, snp)
            if key not in sample_cache:
                sample_cache[key] = next(iter(parse_fasta(sample_path).values()))
            read_seqs = load_reads_by_name(reads_path)
            acc = restoration_accuracy(ref_cache[size], sample_cache[key], result_path, read_seqs)
            restore = f"{acc*100:.2f}%"

        # ── 시간 ──
        elapsed = "N/A"
        time_path = os.path.join(results_dir, subdir, "time.txt")
        if os.path.exists(time_path):
            with open(time_path) as f:
                elapsed = f"{float(f.read().strip()):.1f}s"

        snp_num = int(snp.split("_")[1])
        entries.append((SIZE_BP[size], snp_num, cov, tag, size, f"{snp_num}%", f"{cov}x",
                        f"{mapped/total*100:.1f}%", prec, restore, elapsed))

    entries.sort(key=lambda x: (x[0], x[1], x[2]))

    print(f"\n{'='*92}")
    print(f"{'tag':>8} {'N':>6} {'SNP':>5} {'Cov':>6} "
          f"{'AlignRate':>10} {'Precision':>10} {'RestoreAcc':>11} {'Time':>8}")
    print("-" * 92)

    rows = []
    for _, _, _, t, size, snp_pct, cov_lbl, align, prec, restore, elapsed in entries:
        print(f"{t:>8} {size:>6} {snp_pct:>5} {cov_lbl:>6} "
              f"{align:>10} {prec:>10} {restore:>11} {elapsed:>8}")
        rows.append((t, size, snp_pct, cov_lbl, align, prec, restore, elapsed))

    summary_path = os.path.join(results_dir, "summary.tsv")
    with open(summary_path, "w") as f:
        f.write("tag\tN\tSNP\tCoverage\tAlignRate\tPrecision\tRestoreAcc\tTime\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"  → {summary_path}")


dirs = sys.argv[1:] or [d for d in sorted(os.listdir(".")) if d.startswith("results") and os.path.isdir(d)]
for d in dirs:
    analyze(d.rstrip("/"))
