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
fCO2, fH2O, fCH4, fC3H8 = -394, -286, -75, -105
dP = 3*fCO2 + 4*fH2O - fC3H8; dM = fCO2 + 2*fH2O - fCH4
print(dP, dM, dP + 4*44); assert dP == -2221 and dM == -891 and dP + 4*44 == -2045
for e in ['C3H8 + 5O2 -> 3CO2 + 4H2O', 'CH4 + 2O2 -> CO2 + 2H2O']: assert check_equation(e) is None
x = Symbol('x'); xs = solve(Eq(-dM*x - dP*(1-x), 1423), x)[0]; print('x', xs); assert xs == Rational(3,5)
O2 = xs*2 + (1-xs)*5; air = O2*5*Rational(224,10); print('air L', float(air)); assert rnd(air,3) == 358
