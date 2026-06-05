import time
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.eunho.utils import parse_fasta, stream_reads, write_tsv
from src.eunho.index_builder import build_index
from src.eunho.mapper import map_reads

D = 3  # mismatch 수 


def run_mapping(ref_path, reads_path, output_path, max_repeat=500, d=None):
    d = d if d is not None else D

    # 1 Reference genome 로드
    print(f"[1] Reference genome 로드 from {ref_path}")
    t0 = time.time()
    genome_dict = parse_fasta(ref_path)
    chrom = list(genome_dict.keys())[0]
    genome = genome_dict[chrom]
    print(f"  염기서열: {chrom}, 길이: {len(genome):,} bp  ({time.time()-t0:.2f}s)")

    # 2. k-mer 인덱스 빌드 (Rabin-Karp 롤링 해시을 이용하여 구한다.)
    first_read = next(stream_reads(reads_path), None)
    if first_read is None:
        print("error: reads 파일 없음")
        sys.exit(1)
    L = len(first_read[1])
    k = L // (d + 1) # seed의 길이는 D를 기반으로 설정한다.
    print(f"seed 길이: L={L}, D={d} → k={k}") 

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
