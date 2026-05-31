import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.naive.index_builder import hash_kmer
from src.naive.utils import count_mismatches


def map_read(read, genome, index, k, max_mismatches):
    """최적화 없는 read mapping.

    제거된 최적화:
        - seen set 제거: 같은 후보 위치를 중복 검증
        - 조기 종료 제거: mismatch가 D를 넘어도 끝까지 비교
    """
    L = len(read)
    D = max_mismatches
    n_genome = len(genome)
    hits = []

    # seen set 없음: 중복 검증 허용
    for s in range(D + 1):
        seg_start = s * k
        seg_end   = seg_start + k
        if seg_end > L:
            break
        seed = read[seg_start:seg_end]

        seed_hash  = hash_kmer(seed)
        candidates = index.get(seed_hash, [])

        for ref_seed_pos in candidates:
            read_start = ref_seed_pos - seg_start

            if read_start < 0 or read_start + L > n_genome:
                continue

            # seen 체크 없음

            ref_sub = genome[read_start:read_start + L]
            mm = count_mismatches(ref_sub, read)   # 조기 종료 없음

            if mm <= D:
                hits.append((read_start, mm))

    hits.sort(key=lambda x: x[1])
    return hits


def map_reads(reads, genome, index, k, max_mismatches):
    for read_name, read_seq in reads:
        hits = map_read(read_seq, genome, index, k, max_mismatches)
        yield read_name, hits
