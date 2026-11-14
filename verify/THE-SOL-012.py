# 検証スクリプト（g13 作成）：python3 -I このファイル（作業ディレクトリのモジュールに依存しない）
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, Symbol, solve, Eq, nsimplify
ID = __file__.replace('\\', '/').split('/')[-1][:-3]
ITEM = json.load(open(f'/home/claude/chem/out/items/{ID}/item.json', encoding='utf-8'))
ARR = re.compile(r'<=>|->')


def ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0:
            return out
        k, d = j + 4, 1
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0)
            k += 1
        out.append(s[j + 4:k - 1])
        i = k


def check_all_equations():
    """問題文・答え・解説の \\ce{…} のうち矢印1つの完全な式をすべて検査"""
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(ITEM[key]):
            if len(ARR.findall(b)) != 1 or '\\boxed' in b or re.search(r'(?<![A-Za-z])[Rn](?![a-z])', b):
                continue
            r = check_equation(b)
            assert r is None, (key, b, r)
            n += 1
    print(f'反応式 {n} 本：原子数・電荷のつり合い OK')


def sig(x, n):
    """有効数字 n 桁に丸めた値"""
    return float(f'{float(x):.{n}g}')


def has(key, *frags):
    for f in frags:
        assert f in ITEM[key], (key, f)


def n_answers():
    a = re.findall(r'(?:^|(?<=[\s　]))\((\d+)\)', ITEM['answer'])
    k = ITEM['problem'].count('\\item')
    assert a == [str(i) for i in range(1, k + 1)], (a, k)
    print(f'小問 {k} 個と答えの番号が一致')

# ---------------------------------------------------------------- 本体
check_all_equations()
n_answers()
K = {20: 32, 40: 64, 60: 109, 80: 169}; N = {20: 36, 40: 36, 60: 37, 80: 38}
# (1) 90 g/100 g となる温度
lo = max(t for t in K if K[t] < 90); hi = min(t for t in K if K[t] >= 90)
assert (lo, hi) == (40, 60) and all(N[t] > 10 for t in N)
# (2)
p2 = 90 - K[20]; assert p2 == 58 and 10 < N[20]
# (3)
W = Q(90) * 100 / K[80]; assert W > Q(10) * 100 / N[80]
print('Wmin', float(W)); assert round(float(W)) == 53
p3 = 90 - W * K[20] / 100; print('ppt3', float(p3)); assert round(float(p3)) == 73
assert W * N[20] / 100 > 10
# 53 g（丸めた水量）で計算しても同じ整数
assert round(90 - 53 * 0.32) == 73
# (4)
assert 200 * N[20] / 100 > 50 and 200 * K[20] / 100 > 10      # 冷却では析出なし
p4 = 50 - 50 * N[80] / 100; assert p4 == 31 and 50 * K[80] / 100 > 10
# 実測値（化学便覧）との比較：KNO3 31.6/63.9/109/169，NaCl 35.8/36.3/37.1/38.0
for t, v in {20: 31.6, 40: 63.9, 60: 109.2, 80: 168.8}.items(): assert abs(K[t] - v) / v < 0.02
for t, v in {20: 35.8, 40: 36.3, 60: 37.1, 80: 38.0}.items(): assert abs(N[t] - v) < 0.8
has('answer', '$58\\,\\mathrm g$', '水 $53\\,\\mathrm g$', '$73\\,\\mathrm g$', '$31\\,\\mathrm g$')
print('THE-SOL-012 OK')
