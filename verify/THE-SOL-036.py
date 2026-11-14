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
PA, PB = Rational(136000), Rational(54000)
z = Rational(2, 5); RT = Rational(83, 10)*1000*363
Pb = z*PA + (1 - z)*PB; yb = z*PA/Pb; print('(1)', float(Pb), float(yb)); assert rnd(Pb, 3) == 8.68e4 and rnd(yb, 3) == 0.627
Pd = 1/(z/PA + (1 - z)/PB); xd = z*Pd/PA; print('(2)', float(Pd), float(xd)); assert rnd(Pd, 3) == 7.12e4 and rnd(xd, 3) == 0.209
P3 = Rational(80000); x3 = (P3 - PB)/(PA - PB); y3 = x3*PA/P3
V = (z - x3)/(y3 - x3); L = 1 - V; vol = V*RT/P3
print('(3)', float(x3), float(y3), float(L), float(V), float(vol))
assert rnd(x3, 3) == 0.317 and rnd(y3, 3) == 0.539 and rnd(V, 3) == 0.374 and rnd(L, 3) == 0.626 and rnd(vol, 3) == 14.1
assert Pd < P3 < Pb
# 体積から状態を求める一般の解法（数値）で (3) を逆に確認：体積 vol で圧力 8.00e4
import math
def state_from_volume(Vtot):
    Pall = RT/Vtot
    if Pall <= Pd:
        return float(Pall), 0.0
    lo, hi = float(Pd), float(Pb)
    for _ in range(200):
        P = (lo + hi)/2
        x = (P - PB)/(PA - PB); y = x*PA/P; v = (z - x)/(y - x)
        if float(v*RT/P) > Vtot:
            lo = P
        else:
            hi = P
    return P, float(1 - v)
P_, L_ = state_from_volume(float(vol)); print('vol→P', P_, L_); assert abs(P_ - 8e4) < 1
# (4)
P4, L4 = state_from_volume(50.0); print('(4)', P4, L4); assert L4 == 0 and rnd(P4, 2) == 6.0e4
print('ALL OK')
