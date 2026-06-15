import time
from collections import defaultdict


def parse_fasta(path):
    sequences = {}
    current_name = None
    current_seqs = []

    with open(path) as f:
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


def stream_reads(path):
    with open(path) as f:
        for idx, line in enumerate(f):
            seq = line.rstrip().upper()
            if seq:
                yield f"read_{idx:06d}", seq


def count_mismatches(s1, s2, max_d):
    mm = 0
    for a, b in zip(s1, s2):
        if a != b:
            mm += 1
            if mm > max_d:
                return mm
    return mm


def write_tsv(hits_list, path):
    with open(path, "w") as f:
        f.write("read_name\tchrom\tstart\tmismatches\n")
        for item in hits_list:
            f.write("\t".join(str(x) for x in item) + "\n")


BASE = 131  # 소수 base
MOD = 10**9 + 7  # 소수 mod


def build_index(genome, k, max_repeat=500):
    """Rabin-Karp 롤링 해시로 genome의 k-mer 인덱스를 O(N)에 빌드.

    python 딕셔너리 값으로 반환합니다. hash -> position
    """
    n = len(genome)
    index = defaultdict(list)

    # base^(k-1) mod MOD 미리 계산하여 사용합니다.
    base_pow = pow(BASE, k - 1, MOD)

    # k-mer 해시를 시작하는 부분
    h = 0
    for c in genome[:k]:
        h = (h * BASE + ord(c)) % MOD
    index[h].append(0)

    # 롤링 해시: 이전 해시값을 이용하여 다음 해시 값을 O(1) 시간에 구한다
    for i in range(1, n - k + 1):
        h = (h - ord(genome[i - 1]) * base_pow % MOD + MOD) % MOD  # 이전 해시에서 맨 앞 문자 제거
        # MOD 연산시 음수가 나올 수 있기에 + MOD 후 % MOD를 이용해서 보정하도록 하였다
        h = (h * BASE + ord(genome[i + k - 1])) % MOD  # 자릿 수 보정 뒤 맨 뒤 문자 추가
        index[h].append(i)  # 해시에 저장

    # 반복 서열 필터링: 너무 많이 등장하는 k-mer 제거
    # 이전 실험 조건인 chr1 에선 유효했지만 random하게 만들어진 reference genome의 경우 유의미하지 않음
    if max_repeat > 0:
        to_delete = [key for key, positions in index.items() if len(positions) > max_repeat]
        for key in to_delete:
            del index[key]

    return dict(index)


def hash_kmer(kmer):
    """탐색을 위해 단일 k-mer를 해시하기 위한 메서드"""
    h = 0
    for c in kmer:
        h = (h * BASE + ord(c)) % MOD
    return h


def map_read(read, genome, index, k, max_mismatches):
    L = len(read)
    D = max_mismatches
    n_genome = len(genome)
    hits = []
    seen = set()  # 같은 후보 위치 중복 검증 방지

    # Pigeonhole: (D+1)개 seed 각각을 인덱스에서 조회
    for s in range(D + 1):
        seg_start = s * k
        seg_end = seg_start + k
        seed = read[seg_start:seg_end]

        seed_hash = hash_kmer(seed)
        candidates = index.get(seed_hash, [])

        for ref_seed_pos in candidates:
            # read가 시작될 reference 위치
            # seed를 이용해 read의 시작위치를 계산하여 저장한다.
            read_start = ref_seed_pos - seg_start

            if read_start < 0 or read_start + L > n_genome:
                continue
            # 이미 후보에 있는 read는 건너뛰도록 하여 실행시간을 단축하였다.
            if read_start in seen:
                continue
            seen.add(read_start)

            # 전체 read 검증 trivial 문자 비교를 사용하였다.
            ref_sub = genome[read_start:read_start + L]
            mm = count_mismatches(ref_sub, read, D)

            if mm <= D:
                hits.append((read_start, mm))

    # mismatch 수 기준 정렬
    # 가장 mismatch가 적은 reads를 매핑한다.
    hits.sort(key=lambda x: x[1])
    return hits


# 하나의 reads를 매핑하기 위한 매서드
# reads 하나를 매핑하고 제거하여 메모리 절약
def map_reads(reads_iter, genome, index, k, max_mismatches):
    for read_name, read_seq in reads_iter:
        hits = map_read(read_seq, genome, index, k, max_mismatches)
        yield read_name, hits


def run_mapping(ref_path, reads_path, output_path, max_repeat=500, d=3):
    # 1. Reference genome 로드
    print(f"[1] Reference genome 로드 from {ref_path}")
    t0 = time.time()
    genome_dict = parse_fasta(ref_path)
    chrom = list(genome_dict.keys())[0]
    genome = genome_dict[chrom]
    print(f"  염기서열: {chrom}, 길이: {len(genome):,} bp  ({time.time()-t0:.2f}s)")

    # 2. k-mer 인덱스 빌드 (Rabin-Karp 롤링 해시을 이용하여 구한다.)
    first_read = next(stream_reads(reads_path))
    L = len(first_read[1])
    k = L // (d + 1)  # seed의 길이는 D를 기반으로 설정한다.
    print(f"seed 길이: L={L}, D={d} -> k={k}")

    print("[2] 인덱스 빌드")
    t1 = time.time()
    index = build_index(genome, k, max_repeat=max_repeat)
    t2 = time.time()
    print(f"인덱스 완료: {len(index):,} 고유 k-mer : ({t2-t1:.2f}s)")

    # 3. Read 매핑 (Pigeonhole)
    print(f"[3] Read 매핑 중")
    t3 = time.time()
    results = []
    mapped = 0
    total = 0

    # stream_reads 를 이용해 한번에 하나의 reads 씩 수행.
    # 전체 reads를 메모리에 올리지 않아도 되기에 메모리 절약이 됨.
    for read_name, hits in map_reads(stream_reads(reads_path), genome, index, k, d):
        total += 1
        if hits:
            mapped += 1
            best_pos, best_mm = hits[0]
            results.append((read_name, chrom, best_pos, best_mm))
        else:
            results.append((read_name, "*", -1, -1))

        if total % 10000 == 0:
            elapsed = time.time() - t3
            print(f"{total:,} reads 처리 중 ({total/elapsed:.0f} reads/s)", end="\r")

    elapsed_map = time.time() - t3
    print(f"\n총 {total:,} reads, {mapped:,} 매핑 성공 ({mapped/total*100:.1f}%)  ({elapsed_map:.2f}s)")
    print(f"처리 속도: {total/elapsed_map:,.0f} reads/s")

    # 4. 결과 저장
    print(f"[4] 결과 저장: {output_path}")
    write_tsv(results, output_path)
