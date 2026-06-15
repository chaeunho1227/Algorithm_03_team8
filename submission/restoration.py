def load_reads_by_name(reads_path):
    reads = {}
    with open(reads_path) as f:
        for idx, line in enumerate(f):
            seq = line.rstrip().upper()
            if seq:
                reads[f"read_{idx:06d}"] = seq
    return reads


def restoration_and_snp_recall(reference, sample, result_path, reads):
    N = len(reference)
    coverage = {}

    with open(result_path) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 4 or parts[1] == "*":
                continue
            seq = reads.get(parts[0])
            if seq is None:
                continue
            pos = int(parts[2])
            if pos < 0 or pos + len(seq) > N:
                continue
            for j, base in enumerate(seq):
                coverage[pos + j] = base

    restored = list(reference)
    for i, base in coverage.items():
        restored[i] = base

    N = min(len(restored), len(sample))
    match = truth_snp = detected_snp = 0
    for i in range(N):
        s = sample[i]
        if restored[i] == s:
            match += 1
        if reference[i] != s:
            truth_snp += 1
            if restored[i] == s:
                detected_snp += 1

    acc = match / N if N else 0.0
    recall = detected_snp / truth_snp if truth_snp else 0.0
    return acc, recall
