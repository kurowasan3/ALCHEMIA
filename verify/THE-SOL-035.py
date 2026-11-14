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
nNaCl = Rational(200, 585)            # mol
Lf = Rational(6000, 18)               # J/g
k = Rational(185, 100)*2
m1 = nNaCl/Rational(80, 1000); t1 = -k*m1; print('(1)', float(t1)); assert round(float(t1), 1) == -15.8
# (2) 最終温度 -10.0 → 液相の水
W = nNaCl/(Rational(10)/k)*1000; melt = W - 80; conc = 20/(20 + W)*100
print('(2)', float(melt), float(conc)); assert rnd(melt, 3) == 46.5 and rnd(conc, 3) == 13.7
# (3) 熱収支（エンタルピーの和が一定）を独立に：各物質の基準状態（0 ℃ の水）からのエンタルピー
t0 = symbols('t0')
Hini = 4*100*t0 + 100*(-Lf + 2*0)          # 溶液 100 g（4.0 J/gK），氷 100 g（0 ℃，−Lf）
# 溶液の熱容量は質量×4.0 とし，溶液の質量 = 20 + 液相の水
Hfin = 4*(20 + W)*(-10) + (100 - melt)*(-Lf + 2*(-10))
sol = solve(Eq(Hini, Hfin), t0)[0]; print('(3)', float(sol)); assert round(float(sol), 1) == 21.4
# (4) 氷が融けきる境目
W4 = 180; T4 = -k*nNaCl/Rational(W4, 1000); print('T4', float(T4))
t4 = solve(Eq(4*100*t0 + 100*(-Lf), 4*(20 + W4)*T4), t0)[0]; print('(4)', float(t4)); assert rnd(t4, 2) == 69
# t0 > t4 で氷が残らないこと：t0=80 で氷が残ると仮定すると melt>100 になる
x, T = symbols('x T')
def final_state(tt):
    # 氷 x g が融けて温度 T（凝固点）で共存
    eqs = [Eq(T, -k*nNaCl/((80 + x)/1000)), Eq(4*100*tt + 100*(-Lf), 4*(100 + x)*T + (100 - x)*(-Lf + 2*T))]
    s = nsolve([e.lhs - e.rhs for e in eqs], [x, T], [40, -8])
    return float(s[0]), float(s[1])
print('t0=21.4 →', final_state(sol)); assert abs(final_state(sol)[1] + 10) < 1e-6
xm, Tm = final_state(80); print('t0=80 →', xm, Tm); assert xm > 100
xm, Tm = final_state(60); print('t0=60 →', xm, Tm); assert xm < 100
print('ALL OK')
