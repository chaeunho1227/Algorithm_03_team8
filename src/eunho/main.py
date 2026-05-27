import time
import sys

from .utils import parse_fasta, stream_reads, write_tsv
from .index_builder import build_index
from .mapper import map_reads

# ── 알고리즘 파라미터 ──────────────────────────────────────
D          = 3    # 허용 mismatch 수
MAX_REPEAT = 500  # 인덱스에서 제외할 최대 반복 횟수
# ──────────────────────────────────────────────────────────


def run_mapping(ref_path, reads_path, output_path):
    # 1. Reference genome 로드
    print(f"[1/4] Reference genome 로드 중: {ref_path}")
    t0 = time.time()
    genome_dict = parse_fasta(ref_path)
    chrom = list(genome_dict.keys())[0]
    genome = genome_dict[chrom]
    print(f"  chromosome: {chrom}, 길이: {len(genome):,} bp  ({time.time()-t0:.2f}s)")

    # 2. k-mer 인덱스 빌드 (Rabin-Karp 롤링 해시)
    first_read = next(stream_reads(reads_path), None)
    if first_read is None:
        print("오류: reads 파일이 비어 있습니다.")
        sys.exit(1)
    L = len(first_read[1])
    k = L // (D + 1)
    print(f"  seed 길이: L={L}, D={D} → k={k}")

    print(f"[2/4] k-mer 인덱스 빌드 중 (k={k}, max_repeat={MAX_REPEAT})")
    t1 = time.time()
    index = build_index(genome, k, max_repeat=MAX_REPEAT)
    print(f"  인덱스 완료: {len(index):,} 고유 k-mer  ({time.time()-t1:.2f}s)")

    # 3. Read 매핑 (Pigeonhole)
    print(f"[3/4] Read 매핑 중 (D={D})")
    t2 = time.time()
    results = []
    mapped = 0
    total = 0

    for read_name, hits in map_reads(stream_reads(reads_path), genome, index, k, D):
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
