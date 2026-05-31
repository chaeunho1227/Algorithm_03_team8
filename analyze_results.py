import os, re, sys

SIZE_BP = {"200K": 200_000, "600K": 600_000, "1800K": 1_800_000}
PATTERN = re.compile(r"^(200K|600K|1800K)_(snp_[135])(?:_d\d+)?_cov(\d+)x$")

def analyze(results_dir):
    data_dir = "data/random" if "random" in results_dir else "data"

    print(f"\n{'='*68}\n  {results_dir}\n{'='*68}")
    print(f"{'Config':<32} {'Align':>7} {'Prec':>7} {'Recall':>7}")
    print("-" * 55)

    for subdir in sorted(os.listdir(results_dir)):
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

        print(f"{subdir:<32} {mapped/total*100:>6.1f}% {prec:>7} {rec:>7}")

dirs = sys.argv[1:] or [d for d in sorted(os.listdir(".")) if d.startswith("results") and os.path.isdir(d)]
for d in dirs:
    analyze(d.rstrip("/"))
