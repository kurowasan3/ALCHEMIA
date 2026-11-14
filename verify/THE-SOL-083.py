# THE-SOL-083 検証：霧状の吸収液による気体の吸収（ヘンリーの法則・液滴の大きさと吸収速度）
import sympy as sp
n0, V, L, H, R, T, P, d, k, c = sp.symbols('n0 V L H R T P d k c', positive=True)
# 物質収支：気相 PV/RT ＋ 液相 (P/H)L = n0
Psol = sp.solve(sp.Eq(P*V/(R*T) + P/H*L, n0), P)[0]
print('P =', sp.simplify(Psol))
nabs = sp.simplify(Psol/H*L); print('absorbed =', nabs)
assert sp.simplify(nabs - n0*L*R*T/(H*V + L*R*T)) == 0
# 液滴の直径 d→d/3（全体積 L 一定）
N = 6*L/(sp.pi*d**3); area1 = sp.pi*d**2
rate1 = k*area1*c; total = N*rate1
f = lambda e: sp.simplify(e.subs(d, d/3)/e)
print('最終吸収量の比', f(nabs), ' 1個あたり速度の比', f(rate1), ' 全体の速度の比', f(total))
assert f(nabs) == 1 and f(rate1) == sp.Rational(1, 9) and f(total) == 3
# 数値：V=8.31 L, T=300 K, R=8.31e3, H=3.00e5 Pa·L/mol, L=1.00 L
vals = {V: 8.31, T: 300, R: 8.31e3, H: 3.00e5}
frac1 = (nabs/n0).subs(vals).subs(L, 1.00); frac2 = (nabs/n0).subs(vals).subs(L, 2.00)
print('吸収割合 L=1:', float(frac1), ' L=2:', float(frac2), ' 比', float(frac2/frac1))
assert abs(float(frac1) - 0.5) < 1e-9 and abs(float(frac2) - 2/3) < 1e-9
print('OK')
