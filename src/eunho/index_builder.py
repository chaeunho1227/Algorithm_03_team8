from collections import defaultdict

BASE = 131 # 소수 base
MOD  = 10**9 + 7  # 소수 mod


def build_index(genome, k, max_repeat=500):
    """Rabin-Karp 롤링 해시로 genome의 k-mer 인덱스를 O(N)에 빌드.
    
    python 딕셔너리 값으로 반환합니다. hash -> position
    """
    n = len(genome)
    if n < k:
        return {}

    index = defaultdict(list)

    # base^(k-1) mod MOD 미리 계산하여 사용합니다.
    base_pow = pow(BASE, k - 1, MOD)

    # k-mer 해시를 시작하는 부분
    h = 0
    for c in genome[:k]:
        h = (h * BASE + ord(c)) % MOD
    index[h].append(0)

    # 롤링 해시 : 이전 해시값을 이용하여 다음 해시 값을 O(1) 시간에 구한다
    for i in range(1, n - k + 1):
        h = (h - ord(genome[i - 1]) * base_pow % MOD + MOD) % MOD # 이전 해시에서 맨 앞 문자 제거  
        # MOD 연산시 음수가 나올 수 있기에 + MOD 후 % MOD를 이용해서 보정하도록 하였다
        h = (h * BASE + ord(genome[i + k - 1])) % MOD # 자릿 수 보정 뒤 맨 뒤 문자 추가
        index[h].append(i) #해시에 저장

    # 반복 서열 필터링: 너무 많이 등장하는 k-mer 제거
    # 이전 실헝 조건인 chr1 에선 유효했지만 random하게 만들어진 reference genome의 경우 유의미하지 않음
    if max_repeat > 0:
        to_delete = [key for key, positions in index.items()
                     if len(positions) > max_repeat]
        for key in to_delete:
            del index[key]

    return dict(index)


def hash_kmer(kmer):
    """탐색을 위해 단일 k-mer를 해시하기 위한 매서드"""
    h = 0
    for c in kmer:
        h = (h * BASE + ord(c)) % MOD
    return h
