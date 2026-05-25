import gzip


def parse_fasta(path):
    """FASTA 파일 파싱 → {name: sequence} dict"""
    opener = gzip.open if path.endswith(".gz") else open
    sequences = {}
    current_name = None
    current_seqs = []

    with opener(path, "rt") as f:
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


def stream_fastq(path):
    """FASTQ 파일을 한 번에 하나씩 (name, seq) yield"""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip().upper()
            f.readline()   # '+' 줄
            f.readline()   # quality 줄
            yield header[1:].rstrip().split()[0], seq


def stream_fasta_reads(path):
    """FASTA 형식 reads 파일을 한 번에 하나씩 (name, seq) yield"""
    opener = gzip.open if path.endswith(".gz") else open
    name = None
    seqs = []
    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seqs).upper()
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
        if name is not None:
            yield name, "".join(seqs).upper()


def stream_plain_reads(path):
    """Plain-text reads 파일을 한 번에 하나씩 (name, seq) yield
    팀 레포의 reads_*.txt 형식 (헤더 없이 서열 한 줄씩)을 지원한다."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for idx, line in enumerate(f):
            seq = line.rstrip().upper()
            if seq:
                yield f"read_{idx:06d}", seq


def stream_reads(path):
    """FASTA / FASTQ / plain-text 자동 판별해서 reads yield"""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        first = f.readline()

    if first.startswith("@"):
        return stream_fastq(path)
    elif first.startswith(">"):
        return stream_fasta_reads(path)
    else:
        return stream_plain_reads(path)


def count_mismatches(s1, s2, max_d):
    """두 문자열의 mismatch 수를 trivial 비교로 계산 (조기 종료 포함)
    max_d 초과 시 max_d+1 반환"""
    mm = 0
    for a, b in zip(s1, s2):
        if a != b:
            mm += 1
            if mm > max_d:
                return mm
    return mm



def write_tsv(hits_list, path):
    """매핑 결과를 TSV로 저장
    hits_list: [(read_name, chrom, pos, mismatches), ...]
    """
    with open(path, "w") as f:
        f.write("read_name\tchrom\tstart\tmismatches\n")
        for item in hits_list:
            f.write("\t".join(str(x) for x in item) + "\n")
