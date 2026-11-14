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
n = Rational(280,100)/56; assert n == Rational(1,20)
q1 = 100*Rational(42,10)*Rational(68,10); q2 = 100*Rational(42,10)*Rational(135,10)
d1 = -q1/n/1000; d2 = -q2/n/1000; dn = d2 - d1
print(float(d1), float(d2), float(dn)); assert rnd(d1,2) == -57 and rnd(d2,3) == -113 and rnd(dn,2) == -56
assert Rational(6,10)*Rational(1,10) > n
T0 = (100*25 + 50*22)/Rational(150); 
for dh in (56, -dn):
    T = T0 + Rational(25,1000)*dh*1000/(150*Rational(42,10)); print(float(T0), float(T)); assert round(float(T),1) == 26.2
