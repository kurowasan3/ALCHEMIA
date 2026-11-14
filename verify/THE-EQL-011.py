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
v = Rational(4,100000)/5/50; NO, O2 = Rational(2,1000)/5, Rational(1,1000)/5
k = v/(NO**2*O2); print(float(v), k); assert v == Rational(16,10**8) and k == 5000
print('O2 rate', float(v/2)); assert v/2 == Rational(8,10**8)
assert (2*NO)**2*(2*O2)/(NO**2*O2) == 8
r = {c: a**2*b for c, (a, b) in dict(a=(1,2), b=(Rational(3,2),Rational(3,2)), c=(2,1), d=(Rational(5,2),Rational(1,2))).items()}
print(r); assert max(r, key=r.get) == 'c' and min(r, key=r.get) == 'a' and r['c']/r['a'] == 2
print('conversion', float(Rational(4,100000)/Rational(2,1000)))
