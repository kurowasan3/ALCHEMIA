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
# ダイヤモンド→C(気) ＝ ダイヤモンド→黒鉛 (−2) ＋ 黒鉛→C(気) (+717)
check_all_equations()
assert check_equation('C60 + 60O2 -> 60CO2') is None
assert 60 - 90 + (12+20) == 2 and 60*3/2 == 90
dia = 717 + (-396 - (-394)); print(dia, Rational(dia,2), Rational(717)/Rational(3,2))
assert dia == 715 and Rational(dia,2) == Rational(715,2) and Rational(717)/Rational(3,2) == 478
fC60 = 60*-394 - (-25960); at = 60*717 - fC60; print(fC60, at, float(Rational(at,90))); assert round(at/90) == 452
E = Symbol('E'); e = solve(Eq(2*226 + 498 - 4*E, -911), E)[0]; print('Si-O', float(e)); assert round(float(e)) == 465
