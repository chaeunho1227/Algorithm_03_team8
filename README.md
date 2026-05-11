# Algorithm Team 8

## 데이터셋 생성

```bash
python generate_datasets.py # window
python3 generate_datasets.py # Mac
```

## 데이터 구조

```
data/
└── {200K, 600K, 1800K}/
    ├── reference_{size}.festa
    └── snp_{1, 3, 5}/
        ├── sample.festa      # 정답 서열
        ├── snp.tsv           # 정답 SNP 위치
        ├── reads_{M_10x}.txt # 알고리즘 입력 (10x coverage)
        └── reads_{M_20x}.txt # 알고리즘 입력 (20x coverage)
```
