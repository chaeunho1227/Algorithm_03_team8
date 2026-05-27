

def load_truth(path):
    """truth TSV → {read_name: ref_start} dict"""
    truth = {}
    with open(path) as f:
        next(f)  # 헤더 스킵
        for line in f:
            parts = line.rstrip().split("\t")
            name, start = parts[0], int(parts[1])
            truth[name] = start
    return truth


def load_result(path):
    """result TSV → {read_name: (chrom, start, mismatches)} dict"""
    result = {}
    with open(path) as f:
        next(f)  # 헤더 스킵
        for line in f:
            parts = line.rstrip().split("\t")
            name = parts[0]
            chrom = parts[1]
            start = int(parts[2])
            mm = int(parts[3])
            result[name] = (chrom, start, mm)
    return result


def evaluate(truth, result, tolerance=0):
    """정확도 평가.

    Args:
        tolerance: 정답으로 인정할 위치 오차 범위 (bp)

    Returns:
        dict with metrics
    """
    total = len(truth)
    mapped = 0
    correct = 0
    wrong = 0
    unmapped = 0

    for name, true_start in truth.items():
        if name not in result:
            unmapped += 1
            continue

        chrom, pred_start, mm = result[name]

        if chrom == "*" or pred_start < 0:
            unmapped += 1
            continue

        mapped += 1

        if abs(pred_start - true_start) <= tolerance:
            correct += 1
        else:
            wrong += 1

    precision = correct / mapped if mapped > 0 else 0.0
    recall    = correct / total  if total  > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "total":     total,
        "mapped":    mapped,
        "unmapped":  unmapped,
        "correct":   correct,
        "wrong":     wrong,
        "precision": precision,
        "recall":    recall,
        "f1":        f1,
    }


