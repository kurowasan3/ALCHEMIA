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
Kf = Rational(186, 100)
# (1)
a = symbols('a', positive=True)
al1 = solve(Eq(Kf*Rational(1, 10)*(1+a), Rational(279, 1000)), a)[0]
Ka = Rational(1, 10)*al1**2/(1-al1)
print('(1) alpha', al1, 'Ka', Ka); assert al1 == Rational(1, 2) and rnd(Ka, 2) == 0.050
# (2) 連立
c, al, w = symbols('c al w', positive=True)
sol = solve([Eq(Kf*c*(1+al), Rational(744, 1000)), Eq(c*al**2/(1-al), Ka)], [c, al], dict=True)
sol = [s for s in sol if 0 < s[al] < 1]
assert len(sol) == 1, sol
s = sol[0]; wkg = Rational(1, 100)/s[c]; ice = 100 - wkg*1000
print('(2) alpha', s[al], 'c', s[c], 'ice', float(ice)); assert s[al] == Rational(1, 3) and rnd(ice, 3) == 66.7
# (3) 電離度一定
c3 = Rational(4, 10)/Rational(3, 2); ice3 = 100 - Rational(1, 100)/c3*1000
print('(3) ice', float(ice3)); assert rnd(ice3, 3) == 62.5 and ice > ice3
# (4)
Kap = Rational(1, 10)*al1**2
print('(4) ratio', Kap/Ka); assert Kap/Ka == Rational(1, 2)
cmin = Ka*(1-Rational(5, 100))/Rational(5, 100)**2
print('cmin', float(cmin)); assert rnd(cmin, 2) == 19
# 濃いほど alpha が小さい（単調性）
f = lambda cc: float(nsolve(cc*a**2/(1-a) - Ka, a, 0.3))
assert f(1) > f(5) > f(19)
print('ALL OK')
