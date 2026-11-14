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
Rg = Rational(83, 10); ln2, ln3 = Rational(693, 1000), Rational(110, 100)
k300 = ln2/20; k1a, k2a = k300*Rational(3, 4), k300/4
print('(1)', float(k1a), float(k2a)); assert rnd(k1a, 3) == 2.60e-2 and rnd(k2a, 3) == 8.66e-3
k320 = ln2/5; k1b = k2b = k320/2
inv = Rational(1, 300) - Rational(1, 320); assert inv == Rational(1, 4800)
Ea1 = Rg*(3*ln2 - ln3)/inv; Ea2 = Rg*(3*ln2)/inv
print('(2)', float(Ea1), float(Ea2)); assert rnd(Ea1/1000, 2) == 39 and rnd(Ea2/1000, 2) == 83
# 実の自然対数でも同じ桁
assert rnd(8.3*math.log(8/3)*4800/1000, 2) == 39 and rnd(8.3*math.log(8)*4800/1000, 2) == 83
# (3)
T = symbols('T', positive=True)
ratio = lambda TT: 3*exp((Ea2 - Ea1)/Rg*(1/TT - Rational(1, 300)))
T3 = solve(Eq((Ea2 - Ea1)/Rg*(1/T - Rational(1, 300)), -2*ln3), T)[0]; print('(3)', float(T3)); assert rnd(T3, 3) == 343
print('  実対数での比 at T3', float(3*exp(Rational(4800)*math.log(3)*(1/T3 - Rational(1, 300)))))
assert abs(float(3*exp(4800*math.log(3)*(1/T3 - Rational(1, 300)))) - 1/3) < 1e-9
# (4)
T4 = solve(Eq(ln3*(1 + 4800*(1/T - Rational(1, 300))), 2*ln2), T)[0]; print('(4)', float(T4)); assert rnd(T4, 3) == 295
T4e = 1/(1/300 + (math.log(4)/math.log(3) - 1)/4800); print('  実対数', T4e); assert rnd(T4e, 3) == 295
# 生成物比が途中で一定なことを数値積分で確認
k1, k2 = float(k1a), float(k2a); A_, B_, C_ = 1.0, 0.0, 0.0; dt = 0.01
for step in range(1, 20001):
    dA = (k1 + k2)*A_*dt; B_ += k1*A_*dt; C_ += k2*A_*dt; A_ -= dA
    if step in (100, 2000, 20000):
        assert abs(B_/C_ - 3) < 1e-9
assert abs(A_ - math.exp(-(k1 + k2)*200)) < 1e-3
print('ALL OK')
