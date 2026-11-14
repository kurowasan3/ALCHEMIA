# 検証スクリプト（g03 作成）：python3 -I このファイル
from sympy import *
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ID = __file__.replace('\\', '/').split('/')[-1][:-3]
ITEM = json.load(open(f'/home/claude/chem/out2/items/{ID}/item.json', encoding='utf-8'))
ARR = re.compile(r'<=>|->')
GENERIC = re.compile(r'(?<![A-Za-z])[RnX](?![a-z])')


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
            if len(ARR.findall(b)) != 1 or '\\boxed' in b or GENERIC.search(b):
                continue
            r = check_equation(b)
            assert r is None, (key, b, r)
            n += 1
    print(f'反応式 {n} 本：原子数・電荷のつり合い OK')


def check_subitems():
    """小問の数と答えの番号の数が一致するか"""
    k = ITEM['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', ITEM['answer'])
    assert [int(x) for x in nums] == list(range(1, k + 1)), (k, nums)
    print(f'小問 {k} 個と答えの番号が一致')


def near(a, b, rel=0.005):
    a, b = float(a), float(b)
    assert abs(a - b) <= rel * abs(b) + 1e-12, (a, b)


def rnd(x, n):
    """有効数字 n 桁に丸めた値"""
    x = float(x)
    return float(f'{x:.{n}g}')

# ---------------------------------------------------------------- 本体
check_all_equations(); check_subitems()
k = Rational(83, 10)*1000*300/Rational(415, 100)       # Pa/mol
assert k == 600000
n0 = Rational(18000)/k; wv = Rational(3600)/k
def forward(a):
    """燃焼前のメタン a mol（酸素 n0-a）から，27 ℃ の圧力・水滴の有無・発熱量 kJ"""
    b = n0 - a
    if b >= 2*a:
        burn = a; others = a + (b - 2*a); water = 2*a
    else:
        burn = b/2; others = b/2 + (a - b/2); water = b
    gas_w = min(water, wv)
    P = (others + gas_w)*k
    Q = burn*890 - gas_w*44
    return P, water > wv, Q
# (1) a を細かく走査して P=1.62e4 かつ水滴ありになる a をすべて探す（独立な確認）
sols = []
N = 3000
for i in range(1, N):
    a = n0*Rational(i, N)
    P, drop, Q = forward(a)
    if P == 16200 and drop:
        sols.append(a)
print('(1)', [float(s) for s in sols]); assert sols == [Rational(45, 10000), Rational(21, 1000)]
# (2) 発熱量は同じ
Qs = [forward(a)[2] for a in sols]; print('(2) Q', [float(q) for q in Qs])
assert Qs[0] == Qs[1] and rnd(Qs[0], 3) == 3.74
# 燃えたメタン・CO2・水も同じ
def products(a):
    b = n0 - a
    burn = a if b >= 2*a else b/2
    return burn, burn, 2*burn
assert products(sols[0]) == products(sols[1])
# 一般に：水滴があるとき，圧力が同じなら発熱量も同じ（走査で確認）
# （水滴がないときは，酸素過剰側で圧力が燃焼前と同じになり圧力から何も決まらない）
seen = {}
for i in range(1, N):
    P, d, Q = forward(n0*Rational(i, N))
    if not d:
        continue
    if P in seen:
        assert seen[P] == Q
    seen[P] = Q
# (3)
assert sols[1]/n0 == Rational(7, 10); left = sols[1] - (n0 - sols[1])/2
print('(3) 2回目', left); assert left == Rational(165, 10000)
# (4) 水滴が生じない範囲
nodrop = [Rational(i, N) for i in range(0, N + 1) if not forward(n0*Rational(i, N))[1]]
lo = [x for x in nodrop if x < Rational(1, 2)]; hi = [x for x in nodrop if x > Rational(1, 2)]
print('(4)', float(max(lo)), float(min(hi))); assert max(lo) == Rational(1, 10) and min(hi) == Rational(4, 5)
print('ALL OK')
