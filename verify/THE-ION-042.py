# THE-ION-042 検証：各段階の溶液について電荷均衡＋平衡を厳密に解き，近似値と比較
import sys, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
for eq in ['HCO3- + H+ -> H2CO3', 'Mg^2+ + 2OH- -> Mg(OH)2', 'H+ + OH- -> H2O']:
    assert check_equation(eq) is None, eq
K1 = 4.0e-7; Kw = 1.0e-14; Ksp = 1.8e-11
h = sp.symbols('h', positive=True)
def solve_pH(strong_acid_minus_base, CT):
    # 電荷: h + (陽イオン-陰イオンの差) = [HCO3-] + [OH-] ； strong = [Cl-]-[Na+]... を外から与える
    f = h - strong_acid_minus_base - K1*CT/(h+K1) - Kw/h
    sols = [s for s in sp.nsolve(f, h, 1e-3, prec=30) if True] if False else None
    return sp.nsolve(f, h, 1e-4 if strong_acid_minus_base <= 0 else strong_acid_minus_base, prec=30)
# (1) HCl 4.0e-3 + H2CO3 1.0e-3
h1 = solve_pH(4.0e-3, 1.0e-3); pH1 = -math.log10(float(h1)); print('(1)', pH1)
assert abs(pH1 - 2.40) < 0.005
assert abs(K1*1e-3/4e-3 - 1.0e-7) < 1e-12
# (2) H2CO3 1.0e-3 のみ
h2 = sp.nsolve(h - K1*1e-3/(h+K1) - Kw/h, h, 2e-5, prec=30); pH2 = -math.log10(float(h2)); print('(2)', pH2)
assert abs(pH2 - 4.70) < 0.01
# (3) Na+ 3e-3, Mg 1e-3, OH- 5e-3（沈殿前）→ 沈殿後：電荷 [Na+]+2[Mg]=[OH-]、[Mg][OH]^2=Ksp
Q = 1e-3*(5e-3)**2; assert Q > Ksp
x = sp.symbols('x', positive=True)  # [OH-]
sol = sp.nsolve(3e-3 + 2*Ksp/x**2 - x, x, 3e-3, prec=30)
OH = float(sol); Mg = Ksp/OH**2; pH3 = 14 + math.log10(OH)
print('(3) OH', OH, 'Mg', Mg, 'pH', pH3)
assert abs(pH3 - 11.48) < 0.005 and abs(Mg - 2.0e-6) < 0.05e-6
# (4) Na+ 3e-3 流出，Mg 由来 H+ 2e-3；Cl- 4e-3，炭酸 1e-3：強酸分 = Cl - Na = 1e-3
h4 = solve_pH(1.0e-3, 1.0e-3); pH4 = -math.log10(float(h4)); print('(4)', pH4)
assert abs(pH4 - 3.00) < 0.005
print('OK')
