def parse_fasta(path):
    sequences = {}
    current_name = None
    current_seqs = []

    with open(path, "rt") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seqs).upper()
                current_name = line[1:].split()[0]
                current_seqs = []
            else:
                current_seqs.append(line)
        if current_name is not None:
            sequences[current_name] = "".join(current_seqs).upper()

    return sequences


def load_reads(path):
    """최적화 없는 reads 로드.

    제거된 최적화:
        - stream 방식 제거: 파일 전체를 메모리에 한 번에 로드
    """
    reads = []
    with open(path, "rt") as f:
        for idx, line in enumerate(f):
            seq = line.rstrip().upper()
            if seq:
                reads.append((f"read_{idx:06d}", seq))
    return reads


def count_mismatches(s1, s2):
    """최적화 없는 mismatch 계산.

    제거된 최적화:
        - 조기 종료 제거: mismatch 수와 무관하게 끝까지 비교
    """
    return sum(a != b for a, b in zip(s1, s2))


def write_tsv(hits_list, path):
    with open(path, "w") as f:
        f.write("read_name\tchrom\tstart\tmismatches\n")
        for item in hits_list:
            f.write("\t".join(str(x) for x in item) + "\n")
