# THE-RDX-018 検証：Sn2+ による Hg(II) の逐次還元（HgCl2 → Hg2Cl2 → Hg）と沈殿質量の場合分け
import json, re, sys
from sympy import Rational as Q, symbols, solve, N
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
Hg, Cl = Q(201), Q(355, 10)
Hg2Cl2 = 2 * Hg + 2 * Cl; assert Hg2Cl2 == 473
def mass(n, s):
    """HgCl2 n mol に SnCl2 s mol：沈殿 (Hg2Cl2 mol, Hg mol, 質量)"""
    if s <= n / 2: a, b = s, 0
    elif s <= n: a, b = n - s, 2 * (s - n / 2)
    else: a, b = 0, n
    return a, b, a * Hg2Cl2 + b * Hg
n = Q(1, 100)
print([(float(s), float(mass(n, s)[2])) for s in (0, Q(5, 1000), Q(10, 1000), Q(15, 1000))])
assert mass(n, Q(5, 1000))[2] == Q(2365, 1000) and mass(n, Q(10, 1000))[2] == Q(201, 100) and mass(n, Q(15, 1000))[2] == Q(201, 100)
s = symbols('s')
s1 = solve(473 * s - Q(215, 100), s)[0]; s2 = solve(Hg2Cl2 * (n - s) + Hg * 2 * (s - n / 2) - Q(215, 100), s)[0]
print('s candidates', N(s1, 4), N(s2, 4)); assert 0 < s1 <= n / 2 and n / 2 <= s2 <= n
assert round(float(s1) * 1e3, 2) == 4.55 and round(float(s2) * 1e3, 2) == 8.03
a, b, m = mass(n, s2); print('Hg2Cl2', N(a, 4), N(a * Hg2Cl2, 4), 'Hg', N(b, 4), N(b * Hg, 4))
assert round(float(a) * 1e3, 2) == 1.97 and round(float(b) * 1e3, 2) == 6.06 and round(float(a * Hg2Cl2), 3) == 0.933 and round(float(b * Hg), 2) == 1.22
# (4) n 未知，s = 0.0060 mol，灰黒色の沈殿 2.37 g
sx = Q(6, 1000); nn = symbols('nn', positive=True)
cand = []
c1 = solve(Hg * nn - Q(237, 100), nn)[0]        # s >= n（Hg のみ）
if sx >= c1: cand.append(('Hg only', c1))
c2 = solve(Hg2Cl2 * (nn - sx) + Hg * (2 * sx - nn) - Q(237, 100), nn)[0]   # n/2 <= s <= n
if c2 / 2 <= sx <= c2: cand.append(('mixed', c2))
c3 = 473 * sx                                    # s <= n/2 なら白色のみで質量一定 2.838 g
print('c1', N(c1, 4), 'c2', N(c2, 4), 'white mass', N(c3, 4), cand)
assert len(cand) == 1 and cand[0][0] == 'mixed' and round(float(c2) * 1e2, 2) == 1.03 and c3 != Q(237, 100)
it = json.load(open('/home/claude/chem/out2/items/THE-RDX-018/item.json')); k = 0
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') == 1 and 'boxed' not in eq:
            r = check_equation(eq); assert r is None, (eq, r); k += 1; print('eq ok', eq)
assert k >= 3
print('OK')
