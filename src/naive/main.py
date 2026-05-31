import time
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.naive.utils import parse_fasta, load_reads, write_tsv
from src.naive.index_builder import build_index
from src.naive.mapper import map_reads

D          = 3
MAX_REPEAT = 0


def run_mapping(ref_path, reads_path, output_path):
    # 1. Reference genome 로드
    print(f"[1/4] Reference genome 로드 중: {ref_path}")
    t0 = time.time()
    genome_dict = parse_fasta(ref_path)
    chrom = list(genome_dict.keys())[0]
    genome = genome_dict[chrom]
    print(f"  chromosome: {chrom}, 길이: {len(genome):,} bp  ({time.time()-t0:.2f}s)")

    # 2. k-mer 인덱스 빌드 (최적화 없음)
    reads = load_reads(reads_path)   # 전체 메모리 로드
    if not reads:
        print("오류: reads 파일이 비어 있습니다.")
        sys.exit(1)
    L = len(reads[0][1])
    k = L // (D + 1)
    print(f"  seed 길이: L={L}, D={D} → k={k}")

    print(f"[2/4] k-mer 인덱스 빌드 중 (k={k}, MAX_REPEAT 필터 없음)")
    t1 = time.time()
    index = build_index(genome, k)
    print(f"  인덱스 완료: {len(index):,} 고유 k-mer  ({time.time()-t1:.2f}s)")

    # 3. Read 매핑 (최적화 없음)
    print(f"[3/4] Read 매핑 중 (D={D})")
    t2 = time.time()
    results = []
    mapped = 0
    total = 0

    for read_name, hits in map_reads(reads, genome, index, k, D):
        total += 1
        if hits:
            mapped += 1
            best_pos, best_mm = hits[0]
            results.append((read_name, chrom, best_pos, best_mm))
        else:
            results.append((read_name, "*", -1, -1))

        if total % 10000 == 0:
            elapsed = time.time() - t2
            print(f"  {total:,} reads 처리 중... ({total/elapsed:.0f} reads/s)", end="\r")

    elapsed_map = time.time() - t2
    print(f"\n  총 {total:,} reads, {mapped:,} 매핑 성공 ({mapped/total*100:.1f}%)  ({elapsed_map:.2f}s)")
    print(f"  처리 속도: {total/elapsed_map:,.0f} reads/s")

    # 4. 결과 저장
    print(f"[4/4] 결과 저장: {output_path}")
    write_tsv(results, output_path)

    print("\n=== 완료 ===")
    print(f"  총 reads:   {total:,}")
    print(f"  매핑 성공:  {mapped:,} ({mapped/total*100:.1f}%)")
    print(f"  미매핑:     {total-mapped:,}")
    print(f"  인덱스 빌드: {time.time()-t1-elapsed_map:.2f}s")
    print(f"  매핑:       {elapsed_map:.2f}s")
    print(f"  출력:       {output_path}")
