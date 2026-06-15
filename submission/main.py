import time
from collections import defaultdict

BASE = 131
MOD = 10**9 + 7


def parse_fasta(path):
    sequences = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name:
            sequences[name] = "".join(buf).upper()
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


def build_index(genome, k):
    n = len(genome)
    index = defaultdict(list)
    base_pow = pow(BASE, k - 1, MOD)
    h = 0
    for c in genome[:k]:
        h = (h * BASE + ord(c)) % MOD
    index[h].append(0)
    for i in range(1, n - k + 1):
        h = (h - ord(genome[i - 1]) * base_pow % MOD + MOD) % MOD
        h = (h * BASE + ord(genome[i + k - 1])) % MOD
        index[h].append(i)
    return dict(index)


def hash_kmer(kmer):
    h = 0
    for c in kmer:
        h = (h * BASE + ord(c)) % MOD
    return h


def map_read(read, genome, index, k, D):
    L = len(read)
    n_genome = len(genome)
    hits = []
    seen = set()
    for s in range(D + 1):
        seg_start = s * k
        for ref_pos in index.get(hash_kmer(read[seg_start:seg_start + k]), []):
            start = ref_pos - seg_start
            if start < 0 or start + L > n_genome or start in seen:
                continue
            seen.add(start)
            mm = count_mismatches(genome[start:start + L], read, D)
            if mm <= D:
                hits.append((start, mm))
    hits.sort(key=lambda x: x[1])
    return hits


def run_mapping(ref_path, reads_path, output_path, d=3):
    print(f"[1] Reference 로드: {ref_path}")
    genome_dict = parse_fasta(ref_path)
    chrom = list(genome_dict.keys())[0]
    genome = genome_dict[chrom]
    print(f"  {chrom}: {len(genome):,} bp")

    L = len(next(stream_reads(reads_path))[1])
    k = L // (d + 1)
    print(f"[2] 인덱스 빌드 (k={k})")
    t1 = time.time()
    index = build_index(genome, k)
    print(f"  {len(index):,} 고유 k-mer ({time.time()-t1:.2f}s)")

    print("[3] Read 매핑")
    t3 = time.time()
    results = []
    mapped = total = 0
    for read_name, read_seq in stream_reads(reads_path):
        total += 1
        hits = map_read(read_seq, genome, index, k, d)
        if hits:
            mapped += 1
            results.append((read_name, chrom, hits[0][0], hits[0][1]))
        else:
            results.append((read_name, "*", -1, -1))
        if total % 10000 == 0:
            elapsed = time.time() - t3
            print(f"{total:,} reads ({total/elapsed:.0f}/s)", end="\r")

    elapsed = time.time() - t3
    print(f"\n{total:,} reads, {mapped:,} 매핑 ({mapped/total*100:.1f}%)  {elapsed:.2f}s")

    print(f"[4] 결과 저장: {output_path}")
    write_tsv(results, output_path)
