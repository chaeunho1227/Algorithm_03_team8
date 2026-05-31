from collections import defaultdict

BASE = 131
MOD  = 10**9 + 7


def build_index(genome, k, max_repeat=0):
    """최적화 없는 인덱스 빌드.

    제거된 최적화:
        - 롤링 해시 제거: k-mer마다 O(k)씩 새로 계산 → O(N·k)
        - MAX_REPEAT 필터 없음
    """
    n = len(genome)
    index = defaultdict(list)

    for i in range(n - k + 1):
        kmer = genome[i:i + k]
        h = hash_kmer(kmer)   # 매번 O(k)로 새로 계산
        index[h].append(i)

    # MAX_REPEAT 필터 없음

    return dict(index)


def hash_kmer(kmer):
    h = 0
    for c in kmer:
        h = (h * BASE + ord(c)) % MOD
    return h
