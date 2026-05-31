import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.eunho.utils import count_mismatches
from src.eunho.index_builder import hash_kmer


def map_read(read, genome, index, k, max_mismatches):
    L = len(read)
    D = max_mismatches
    n_genome = len(genome)
    hits = []
    seen = set()   # 같은 후보 위치 중복 검증 방지

    # Pigeonhole: (D+1)개 segment 각각을 인덱스에서 조회
    for s in range(D + 1):
        seg_start = s * k
        seg_end = seg_start + k
        if seg_end > L:
            break
        seed = read[seg_start:seg_end]

        seed_hash = hash_kmer(seed)
        candidates = index.get(seed_hash, [])

        for ref_seed_pos in candidates:
            # read가 시작될 reference 위치
            read_start = ref_seed_pos - seg_start

            if read_start < 0 or read_start + L > n_genome:
                continue
            if read_start in seen:
                continue
            seen.add(read_start)

            # 전체 read 검증: trivial 문자 비교 (조기 종료 포함)
            ref_sub = genome[read_start:read_start + L]
            mm = count_mismatches(ref_sub, read, D)

            if mm <= D:
                hits.append((read_start, mm))

    # mismatch 수 기준 정렬 (최적 hit 우선)
    hits.sort(key=lambda x: x[1])
    return hits


def map_reads(reads_iter, genome, index, k, max_mismatches):
    for read_name, read_seq in reads_iter:
        hits = map_read(read_seq, genome, index, k, max_mismatches)
        yield read_name, hits
