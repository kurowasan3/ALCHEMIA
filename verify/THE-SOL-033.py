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
P0, ps = Rational(320, 100)*1000, Rational(288, 100)*1000
w0, u = Rational(10), Rational(1, 2)
p_init = w0/(w0 + u)*P0; print('(1)', float(p_init)); assert rnd(p_init, 3) == 3.05e3 and p_init > ps


def final(n):
    """y を少しずつ増やす代わりに場合分けを独立に：X が尽きる移動量と蒸気圧が p* になる移動量の小さい方"""
    y = symbols('y')
    ya = solve(Eq((w0 - y)/(w0 + u - y)*P0, ps), y)[0]
    yb = 2*n
    yy = min(ya, yb)
    p = (w0 - yy)/(w0 + u - yy)*P0
    return yy, p, ya


y2, p2, ya = final(Rational(1)); print('(2)', y2, float(p2)); assert y2 == 2 and rnd(p2, 3) == 3.01e3 and p2 > ps
y3, p3, _ = final(Rational(4)); print('(3)', y3, float(p3), 'hyd', y3/2, 'X', 4 - y3/2)
assert y3 == Rational(11, 2) and p3 == ps and y3/2 == Rational(11, 4) and 4 - y3/2 == Rational(5, 4)
print('threshold n', ya/2); assert ya/2 == Rational(11, 4)
# 数値シミュレーション（小刻みに水を移す）でも同じ結果
for n, exp in ((1.0, 3010), (4.0, 2880)):
    y, X = 0.0, n
    while X > 1e-9 and (10 - y)/(10.5 - y)*3200 > 2880 + 1e-6:
        dy = min(1e-4, 2*X); y += dy; X -= dy/2
    p = (10 - y)/(10.5 - y)*3200
    assert abs(p - exp) < 5, (n, p)
# (4)
assert 2500 < ps < 3000 and ps/P0 == Rational(9, 10)
print('ALL OK')
