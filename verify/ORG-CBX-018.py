# ORG-CBX-018 検証：酸価・けん化価・ヨウ素価から遊離脂肪酸の割合と脂肪酸の種類を決める
import json, re, sys
from sympy import Rational as Q, symbols, solve, N
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
H, C, O, K, I = Q(1), Q(12), Q(16), Q(39), Q(127)
KOH = K + O + H; I2 = 2 * I
AV, SV, IV = Q(400, 10), Q(193), Q(175)        # mg/g, mg/g, g/100 g
y = AV / 1000 / KOH                     # 遊離脂肪酸 mol（1.00 g 中）
s = SV / 1000 / KOH                     # 3x + y
x = (s - y) / 3
print('y', N(y, 4), 'x', N(x, 4), '3x+y', N(s, 4))
assert round(float(y) * 1e4, 2) == 7.14 and round(float(x) * 1e4, 2) == 9.11
Ma = symbols('Ma')
ma = solve(x * (3 * Ma + 92 - 54) + y * Ma - 1, Ma)[0]   # トリグリセリド M = 3Ma + 92 - 3*18
print('Ma', N(ma, 5)); assert round(float(ma)) == 280
k_from_M = (284 - 280) / 2                        # C17H(35-2k)COOH : M = 284 - 2k
k_from_IV = IV / 100 / I2 / s
print('k(M)', k_from_M, 'k(IV)', N(k_from_IV, 4)); assert k_from_M == 2 and abs(float(k_from_IV) - 2) < 0.01
MT = 3 * 280 + 92 - 54; assert MT == 878 == 57 * 12 + 98 * 1 + 6 * 16   # C57H98O6
wT = x * MT; print('T mass %', N(100 * wT, 4)); assert round(float(100 * wT), 1) == 80.0
frac = y / s; print('hydrolysed %', N(100 * frac, 4)); assert round(float(100 * frac), 1) == 20.7
nH2 = 2 * s; VH2 = nH2 * Q(224, 10); print('H2 mol', N(nH2, 4), 'L', N(VH2, 4)); assert round(float(VH2), 3) == 0.154
newmass = 1 + nH2 * 2 * H
print('new AV', N(AV / newmass, 4), 'new SV', N(SV / newmass, 4)); assert round(float(AV / newmass), 1) == 39.5 and round(float(SV / newmass)) == 190
# 反応式
it = json.load(open('/home/claude/chem/out2/items/ORG-CBX-018/item.json')); n = 0
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') == 1 and 'boxed' not in eq and 'RCOO' not in eq:
            r = check_equation(eq); assert r is None, (eq, r); n += 1; print('eq ok', eq)
assert n >= 3
print('OK')
