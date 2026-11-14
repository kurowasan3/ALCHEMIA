# 検証スクリプト（g03 作成）：python3 -I このファイル
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import *
import re
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


def near(a, b, rel=0.005):
    a, b = float(a), float(b)
    assert abs(a - b) <= rel * abs(b) + 1e-12, (a, b)


def rnd(x, n):
    """有効数字 n 桁に丸めた値"""
    x = float(x)
    return float(f'{x:.{n}g}')

# ---------------------------------------------------------------- 本体
check_all_equations()
d1 = Rational(69,10)*Rational(5,100)/(Rational(95,100)*Rational(128,1000))
d2 = 8*Rational(5,100)/(Rational(95,100)*Rational(154,1000))
print(float(Rational(805,10)-d1), float(69-d2))
assert round(float(Rational(805,10)-d1),1) == 77.7 and round(float(69-d2),1) == 66.3
s1, s2 = d1*20, d2*20
x = Symbol('x'); xe = solve(Eq(Rational(805,10)-s1*x, 69-s2*(1-x)), x)[0]; Te = Rational(805,10)-s1*xe
print('eutectic', float(xe), float(Te)); assert rnd(xe,2) == 0.59 and rnd(Te,2) == 47
T3 = Rational(805,10)-s1*Rational(3,10); print('start', float(T3)); assert rnd(T3,2) == 63
liq = Rational(3,10)/Rational(59,100); fr = (Rational(7,10)-(liq-Rational(3,10)))/Rational(7,10)
print('frac', float(fr)); assert rnd(fr*100,2) == 70
