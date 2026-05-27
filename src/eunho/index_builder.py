from collections import defaultdict

BASE = 131       # 소수 base (DNA 4종류보다 크게 잡아 충돌 감소)
MOD  = 10**9 + 7  # 소수 mod


def build_index(genome, k, max_repeat=500):
    """Rabin-Karp 롤링 해시로 genome의 k-mer 인덱스를 O(N)에 빌드.

    Args:
        genome: reference 문자열 (uppercase)
        k: seed 길이
        max_repeat: 이 이상 등장하는 k-mer는 인덱스에서 제외 (반복서열 필터)

    Returns:
        dict: {hash_value(int) -> [position, ...]}
    """
    n = len(genome)
    if n < k:
        return {}

    index = defaultdict(list)

    # base^(k-1) mod MOD 미리 계산
    base_pow = pow(BASE, k - 1, MOD)

    # 첫 k-mer 해시
    h = 0
    for c in genome[:k]:
        h = (h * BASE + ord(c)) % MOD
    index[h].append(0)

    # 롤링 해시: O(1)씩 업데이트
    for i in range(1, n - k + 1):
        h = (h - ord(genome[i - 1]) * base_pow % MOD + MOD) % MOD # 이전 해시에서 맨 앞 문자 제거 
        # MOD 연산시 음수가 나올 수 있기에 + MOD 후 % MOD를 이용해서 보정한다.
        h = (h * BASE + ord(genome[i + k - 1])) % MOD # 맨 뒤 문자 추가
        index[h].append(i) #해시에 저장

    # 반복 서열 필터링: 너무 많이 등장하는 k-mer 제거
    if max_repeat > 0:
        to_delete = [key for key, positions in index.items()
                     if len(positions) > max_repeat]
        for key in to_delete:
            del index[key]

    return dict(index)


def hash_kmer(kmer):
    """단일 k-mer 문자열의 Rabin-Karp 해시값 계산"""
    h = 0
    for c in kmer:
        h = (h * BASE + ord(c)) % MOD
    return h
