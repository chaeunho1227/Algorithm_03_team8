# Algorithm Team 8 — DNA 서열 복원 (SNP 탐지)

Reference genome에 reads를 매핑하여 SNP를 검출하는 파이프라인.
Pigeonhole Principle 기반 seed-and-extend 매핑 알고리즘을 구현한다.

실행 코드는 모두 `submission/`에 자체 완결적으로 들어 있으며, **random 서열 데이터셋만** 사용한다.

## 실행

```bash
cd submission

# 1. 랜덤 데이터셋 생성 (reference / sample / reads / truth)
python3 generate_random_datasets.py

# 2. 실험 일괄 실행 (3 genome × 3 SNP × 2 coverage = 18조합 → results/ 생성)
python3 run_experiments.py
```

> 두 스크립트 모두 상대 경로(`data/`)를 사용하므로 반드시 `submission/` 안에서 실행한다.

## 파일 구성

```
submission/
├── generate_random_datasets.py  # 랜덤 reference·sample·reads·truth 생성
├── main.py                      # 매핑 파이프라인 (FASTA 파싱 + 인덱스 빌드 + Pigeonhole 매핑)
├── restoration.py               # consensus 복원 + 복원정확도 / SNP 재현율 평가
└── run_experiments.py           # 18개 조건 일괄 실행
```

## 데이터 구조

```
submission/data/
└── {200K, 600K, 1800K}/
    ├── reference_{size}.festa        # 랜덤 ACGT reference
    └── snp_{1, 3, 5}/
        ├── sample.festa             # SNP 주입한 정답 서열
        ├── reads_{M}.txt            # 매핑 입력 (M = genome_size × coverage / 100)
        └── truth_{M}.tsv            # 각 read의 정답 시작 위치
```

## 측정 지표 (4지표)

정밀도는 random reference에서 항상 100%라 제외한다.

- **정렬률** = 정렬된 read 수 / 전체 read 수
- **SNP 재현율** = 탐지된 SNP 위치 수 / 실제 SNP 위치 수
- **복원정확도** = 일치 염기 수 / 전체 염기 수
- **시간** = 매핑 실행 시간

결과는 `submission/results/{size}_{snp}_{cov}x/result.tsv`와 `summary.tsv`에 저장된다.
