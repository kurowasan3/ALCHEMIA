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
Ka = Rational(64, 10)*Rational(1, 10**5)
assert Ka/Rational(1, 10) == Rational(64, 100000)
KD, K2, c = symbols('KD K2 c', positive=True)
org = lambda c_, KD_, K2_: KD_*c_ + 2*K2_*(KD_*c_)**2
sol = solve([Eq(org(Rational(1, 100), KD, K2), Rational(1, 10)), Eq(org(Rational(2, 100), KD, K2), Rational(36, 100))], [KD, K2], dict=True)
print('(2)', sol); assert len(sol) == 1 and sol[0][KD] == 2 and sol[0][K2] == 100
KD0, K20 = 2, 100
dim1 = 2*K20*(KD0*Rational(1, 100))**2/Rational(1, 10); print('dimer %', dim1*100); assert dim1 == Rational(4, 5)
# (3)
r = [s for s in solve(Eq(Rational(1, 10)*(c + org(c, KD0, K20)), Rational(5, 1000)), c) if s > 0]
assert len(r) == 1; c3 = r[0]; f3 = 1 - Rational(1, 10)*c3/Rational(5, 1000)
print('(3)', c3, float(f3*100)); assert c3 == Rational(1, 160) and rnd(f3*100, 3) == 87.5
r2 = [s for s in solve(Eq(Rational(1, 10)*(c + org(c, KD0, K20)), Rational(10, 1000)), c) if s > 0][0]
f3b = 1 - Rational(1, 10)*r2/Rational(10, 1000); print('  2倍', float(f3b)); assert f3b > f3
# (4)
H = Ka
r4 = [s for s in solve(Eq(Rational(1, 10)*(c*(1 + Ka/H) + org(c, KD0, K20)), Rational(5, 1000)), c) if s > 0][0]
f4 = 1 - Rational(1, 10)*r4*(1 + Ka/H)/Rational(5, 1000); print('(4)', float(r4), float(f4*100)); assert rnd(f4*100, 2) == 77
c4r = 4*Rational(232, 100)/1600; f4r = 1 - Rational(1, 10)*2*c4r/Rational(5, 1000); print('  √11=3.32 使用', float(f4r*100)); assert rnd(f4r*100, 2) == 77
H8 = Rational(1, 10**8)
r8 = [s for s in solve(Eq(Rational(1, 10)*(c*(1 + Ka/H8) + org(c, KD0, K20)), Rational(5, 1000)), c) if s > 0][0]
f8 = 1 - Rational(1, 10)*r8*(1 + Ka/H8)/Rational(5, 1000); print('pH8 org %', float(f8*100)); assert f8 < 0.001
# pH 1 で安息香酸自身の H+ を含めても電離は無視できる
print('ALL OK')
