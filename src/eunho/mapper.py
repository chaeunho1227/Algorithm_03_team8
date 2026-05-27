from .utils import count_mismatches
from .index_builder import hash_kmer


def map_read(read, genome, index, k, max_mismatches):
    """Pigeonhole Principle 기반 read mapping.

    핵심 원리:
        read에 최대 D개 mismatch가 있으면, read를 (D+1)개 segment로 나눴을 때
        적어도 1개 segment는 reference에 완전히 일치해야 한다.

    Args:
        read: 매핑할 read 서열
        genome: reference genome 문자열
        index: build_index()로 만든 k-mer 인덱스
        k: seed(segment) 길이
        max_mismatches: 허용 mismatch 수 D

    Returns:
        list of (ref_start, mismatches)  — mismatch 오름차순 정렬
    """
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
    """reads_iter의 모든 read를 매핑하고 결과를 yield.

    Yields:
        (read_name, hits) where hits = [(ref_start, mismatches), ...]
    """
    for read_name, read_seq in reads_iter:
        hits = map_read(read_seq, genome, index, k, max_mismatches)
        yield read_name, hits
