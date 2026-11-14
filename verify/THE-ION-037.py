# THE-ION-037 検証：AgCl・AgI のアンモニア水への溶解（沈殿平衡と錯形成平衡の連立）
import json, re, sys
from sympy import Rational as Q, symbols, solve, sqrt, N, nsolve
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
KspCl, KspI, beta = Q(18, 10) * Q(10)**-10, Q(80, 10) * Q(10)**-17, Q(16, 10) * Q(10)**7
KCl, KI = KspCl * beta, KspI * beta
print('K_AgCl', N(KCl, 3), 'K_AgI', N(KI, 3)); assert round(float(KCl) * 1e3, 2) == 2.88 and round(float(KI) * 1e9, 2) == 1.28
# (2) 1.00 mol/L NH3 に対する AgCl の溶解度（NH3 の消費を考慮）
s = symbols('s', positive=True)
sol = [r for r in solve(s**2 - KCl * (1 - 2 * s)**2, s) if 0 < r < Q(1, 2)][0]
approx = sqrt(KCl) * 1
print('s exact', N(sol, 4), 'approx', N(approx, 4), 'err %', N(100 * (approx - sol) / sol, 3))
assert round(float(sol), 4) == 0.0485 and round(float(approx), 3) == 0.054
# (3) AgCl 1.0e-3 mol + AgI 1.0e-3 mol，100 mL：AgCl を溶かしきる最小の NH3 全濃度
c = Q(1, 100)          # [Cl-]=[Ag(NH3)2+]=0.010（AgI の溶解を無視）
free = c / sqrt(KCl); tot = free + 2 * c
print('free NH3', N(free, 4), 'total', N(tot, 4)); assert round(float(tot), 2) == 0.21 and round(float(free), 3) == 0.186
Ag = KspCl / c; I = KspI / Ag
print('[Ag+]', N(Ag, 3), '[I-] from AgI', N(I, 3)); assert I < Q(1, 10**8)
# (4) 6.0 mol/L NH3 100 mL に溶ける AgI
sI = sqrt(KI) * 6; print('AgI solubility in 6 M', N(sI, 3), 'mol in 0.1 L', N(sI / 10, 3), 'fraction', N(sI / 10 / Q(1, 1000), 3))
assert round(float(sI) * 1e4, 1) == 2.1 and round(float(sI / 10 / Q(1, 1000)) * 100) == 2
# (5) ろ液（[Ag(NH3)2+]=0.010，[NH3]=0.186）に [I-]=1.0e-3 → AgI 沈殿
Agf = c / (beta * free**2); Qp = Agf * Q(1, 1000)
print('[Ag+] filtrate', N(Agf, 3), 'Q', N(Qp, 3)); assert Qp > KspI and round(float(Agf) * 1e8, 1) == 1.8
# 臭化物：Ksp(AgBr)=5.0e-13 で同じろ液に [Br-]=1.0e-3 → 1.8e-11 > 5.0e-13 で沈殿
it = json.load(open('/home/claude/chem/out2/items/THE-ION-037/item.json')); k = 0
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') + eq.count('<=>') == 1 and 'boxed' not in eq:
            r = check_equation(eq); assert r is None, (eq, r); k += 1; print('eq ok', eq)
assert k >= 3
print('OK')
