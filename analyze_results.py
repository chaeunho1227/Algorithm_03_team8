import os, re, sys

SIZE_BP = {"200K": 200_000, "600K": 600_000, "1800K": 1_800_000}
PATTERN = re.compile(r"^(200K|600K|1800K)_(snp_[135])(?:_d\d+)?_cov(\d+)x$")

def analyze(results_dir):
    data_dir = "data/random" if "random" in results_dir else "data"
    tag = "random" if "random" in results_dir else "chr1"

    entries = []
    for subdir in os.listdir(results_dir):
        m = PATTERN.match(subdir)
        result_path = os.path.join(results_dir, subdir, "result.tsv")
        if not m or not os.path.exists(result_path):
            continue

        size, snp, cov = m.group(1), m.group(2), int(m.group(3))
        reads = SIZE_BP[size] * cov // 100
        truth_path = os.path.join(data_dir, size, snp, f"truth_{reads}.tsv")

        total = mapped = 0
        with open(result_path) as f:
            next(f)
            for line in f:
                total += 1
                if line.split("\t")[1] != "*":
                    mapped += 1

        prec = rec = "N/A"
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
            rec  = f"{correct/len(truth)*100:.1f}%"

        snp_num = int(snp.split("_")[1])
        entries.append((SIZE_BP[size], snp_num, cov, tag, size, f"{snp_num}%", f"{cov}x",
                        f"{mapped/total*100:.1f}%", prec, rec))

    entries.sort(key=lambda x: (x[0], x[1], x[2]))

    print(f"\n{'='*78}")
    print(f"{'tag':>8} {'N':>6} {'SNP':>5} {'Coverage':>8} {'AlignRate':>10} {'Precision':>10} {'Recall':>8}")
    print("-" * 78)

    rows = []
    for _, _, _, t, size, snp_pct, cov_lbl, align, prec, rec in entries:
        print(f"{t:>8} {size:>6} {snp_pct:>5} {cov_lbl:>8} {align:>10} {prec:>10} {rec:>8}")
        rows.append((t, size, snp_pct, cov_lbl, align, prec, rec))

    summary_path = os.path.join(results_dir, "summary.tsv")
    with open(summary_path, "w") as f:
        f.write("tag\tN\tSNP\tCoverage\tAlignRate\tPrecision\tRecall\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"  → {summary_path}")

dirs = sys.argv[1:] or [d for d in sorted(os.listdir(".")) if d.startswith("results") and os.path.isdir(d)]
for d in dirs:
    analyze(d.rstrip("/"))
