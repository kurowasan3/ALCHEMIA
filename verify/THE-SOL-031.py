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
RT = Rational(83, 10)*1000*300
p, a = Rational(1, 1000), Rational(1, 100)
# (1)
P1 = (p + a)*RT; print('(1)', float(P1)); assert rnd(P1, 2) == 2.7e4
# (2) ドナン平衡を独立に解く（未知数4つ＋条件）
NaI, ClI, NaII, ClII = symbols('NaI ClI NaII ClII', positive=True)
eqs = [Eq(NaI, a + ClI),            # I 室の電気的中性（Q の負電荷 10×p = a）
       Eq(NaII, ClII),              # II 室の電気的中性
       Eq(ClI + ClII, a),           # Cl の保存（体積等しい）
       Eq(NaI*ClI, NaII*ClII)]      # ドナン平衡
sol = solve(eqs, [NaI, ClI, NaII, ClII], dict=True)
assert len(sol) == 1, sol
s = sol[0]; print('(2)', {k: float(v) for k, v in s.items()})
assert rnd(s[NaI], 3) == 1.33e-2 and rnd(s[ClI], 3) == 3.33e-3 and rnd(s[NaII], 3) == 6.67e-3
# Na の保存も確認
assert s[NaI] + s[NaII] == 2*a
# (3)
dc = (p + s[NaI] + s[ClI]) - (s[NaII] + s[ClII]); P3 = dc*RT
print('(3)', float(P3)); assert rnd(P3, 2) == 1.1e4
# (4)
Mapp = 50/(P1/RT); print('(4) Mapp', float(Mapp)); assert rnd(Mapp, 2) == 4.5e3
c0, x = symbols('c0 x', positive=True)
xs = solve(Eq((a + x)*x, (c0 - x)**2), x)[0]
dcg = simplify((p + a + 2*xs) - 2*(c0 - xs))
print('dc general', dcg)
assert simplify(dcg - (p + a**2/(a + 2*c0))) == 0
assert simplify(dcg.subs(c0, 0) - (p + a)) == 0 and simplify(dcg.subs(c0, a) - dc) == 0
lim = limit(dcg, c0, oo)*RT; print('limit', float(lim)); assert rnd(lim, 2) == 2.5e3
assert rnd(50/p, 2) == 5.0e4
print('ALL OK')
