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
M = 95 + 6*18
assert M == 203
n = Rational(203, 10) / M
w = Rational(1892, 10) + 6*n*18
assert w == 200
m = n / (w/1000)
a = solve(Eq(Rational(240,100), Rational(185,100)*m*(1+2*Symbol('a'))), Symbol('a'))[0]
print('alpha', float(a)); assert rnd(a, 2) == 0.80
m2 = n / Rational(1892, 10000)
a2 = solve(Eq(Rational(240,100), Rational(185,100)*m2*(1+2*Symbol('a'))), Symbol('a'))[0]
print('alpha wrong', float(a2)); assert rnd(a2, 2) == 0.73 and a2 < a
wl = w * Rational(240, 480)
print('ice', float(w - wl)); assert w - wl == 100
print('MgCl2 g/100g water at -4.8C', float(n*95/wl*100))
