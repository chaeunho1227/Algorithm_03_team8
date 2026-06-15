"""복원 정확도 (Restoration accuracy) 계산.

정의: 일치 염기 수 / 전체 염기 수

매핑된 read들로 reference 좌표 위에 복원 게놈을 구성한 뒤,
정답 sample 게놈과 base 단위로 비교한다.

- 커버된 위치: 해당 위치를 덮는 read의 염기로 교체 (마지막 read 기준)
- 커버되지 않은 위치: reference 염기 그대로 둔다 (복원 불가 → 변이 미검출)

복원정확도 = (복원 게놈[i] == sample[i] 인 위치 수) / N
"""


def load_reads_by_name(reads_path):
    reads = {}
    with open(reads_path) as f:
        for idx, line in enumerate(f):
            seq = line.rstrip().upper()
            if seq:
                reads[f"read_{idx:06d}"] = seq
    return reads


def restore_genome(reference, result_path, reads):
    """매핑 결과 + reads 로 복원 게놈을 구성.

    Returns: 복원된 게놈 문자열 (reference 와 동일 길이).
    """
    N = len(reference)
    coverage = {}  # {위치: 염기} — 커버된 위치만 보관

    with open(result_path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            name, chrom, start = parts[0], parts[1], parts[2]
            if chrom == "*":
                continue
            seq = reads.get(name)
            if seq is None:
                continue
            try:
                pos = int(start)
            except ValueError:
                continue
            if pos < 0 or pos + len(seq) > N:
                continue
            for j, base in enumerate(seq):
                coverage[pos + j] = base

    restored = list(reference)
    for i, base in coverage.items():
        restored[i] = base
    return "".join(restored)


def restoration_and_snp_recall(reference, sample, result_path, reads):
    """복원정확도와 SNP 재현율을 한 번의 복원으로 동시 계산.

    - 복원정확도 = 일치 염기 수 / 전체 염기 수
    - SNP 재현율  = 탐지된 SNP 위치 수 / 실제 SNP 위치 수
        실제 SNP 위치 = reference[i] != sample[i] 인 위치
        탐지된 SNP    = 복원 게놈이 그 위치에서 sample 염기로 올바르게 복원된 위치
                        (restored[i] == sample[i] 이고 reference[i] != sample[i])

    Returns: (restoration_accuracy, snp_recall)
    """
    restored = restore_genome(reference, result_path, reads)
    N = min(len(restored), len(sample))

    match = 0
    truth_snp = 0
    detected_snp = 0
    for i in range(N):
        s = sample[i]
        if restored[i] == s:
            match += 1
        if reference[i] != s:            # 실제 SNP 위치
            truth_snp += 1
            if restored[i] == s:         # 올바르게 복원(탐지)된 SNP
                detected_snp += 1

    acc = match / N if N else 0.0
    recall = detected_snp / truth_snp if truth_snp else 0.0
    return acc, recall
