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
NO, O2, N2O2, k1, km1, k2 = symbols('NO O2 N2O2 k1 km1 k2', positive=True)
# 定常状態近似（一般式）から極限を確認
ss = solve(Eq(k1*NO**2, km1*N2O2 + k2*N2O2*O2), N2O2)[0]
v = simplify(2*k2*ss*O2); print('一般の v', v)
assert simplify(v - 2*k1*k2*NO**2*O2/(km1 + k2*O2)) == 0
K = k1/km1
assert simplify(limit(v*km1/k1, km1, oo) - 2*k2*NO**2*O2) == 0     # k_-1 大 → 2 k2 K [NO]^2[O2]
assert simplify(limit(v, O2, oo) - 2*k1*NO**2) == 0                  # [O2] 大 → 2 k1 [NO]^2
# (2)
Ea = -40 + 30; print('(2) Ea', Ea); assert Ea == -10
fac = math.exp(-Ea*1000/8.3*(1/400 - 1/300)); print('  倍率', fac); assert rnd(fac, 2) == 0.37
expo = 1.0e4/8.3/1200; assert abs(expo - 1.0) < 0.01
# (4) 機構 II の一般式と極限
NO3, k1p, km1p, kp = symbols('NO3 k1p km1p kp', positive=True)
ss2 = solve(Eq(k1p*NO*O2, km1p*NO3 + kp*NO3*NO), NO3)[0]
v2 = simplify(2*kp*ss2*NO)
assert simplify(limit(v2*km1p/k1p, km1p, oo) - 2*kp*NO**2*O2) == 0
assert simplify(limit(v2/(NO*O2), NO, oo) - 2*k1p) == 0
# 機構 II では [O2] を大きくしても O2 について1次のまま（O2 で飽和しない）
assert simplify(diff(v2, O2) - v2/O2) == 0
# 機構 I では [NO] を大きくしても NO について2次のまま
assert simplify(diff(v, NO) - 2*v/NO) == 0
print('ALL OK')
