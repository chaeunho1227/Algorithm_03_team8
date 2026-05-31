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
    total   = len(truth)
    mapped  = 0
    correct = 0

    for name, true_start in truth.items():
        if name not in result:
            continue
        chrom, pred_start = result[name]
        if chrom == "*":
            continue
        mapped += 1
        if pred_start == true_start:
            correct += 1

    precision = correct / mapped if mapped > 0 else 0.0
    recall    = correct / total  if total  > 0 else 0.0

    return {"precision": precision, "recall": recall}
