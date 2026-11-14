# THE-THC-046 検証：金属の燃焼の到達状態（相変化を含む断熱の熱収支，窒素による希釈）
import sympy as sp
R = sp.Rational
Q = 250  # kJ (1 mol MO あたり)
Cs, Cl, Cg = R(50), R(60), R(40)   # J/(K mol)
Tm, Tb, T0 = 1200, 2000, 25
Hf, Hv = 40, 250
CN2 = R(30)
def state(nN2):
    a = (Cs + nN2*CN2)*(Tm - T0)/1000
    if Q < a: return ('solid', T0 + Q*1000/(Cs + nN2*CN2))
    if Q < a + Hf: return ('solid+liquid', Tm, (Q - a)/Hf)
    b = a + Hf + (Cl + nN2*CN2)*(Tb - Tm)/1000
    if Q < b: return ('liquid', Tm + (Q - a - Hf)*1000/(Cl + nN2*CN2))
    if Q < b + Hv: return ('liquid+gas', Tb, (Q - b)/Hv)
    return ('gas', None)
s1 = state(0); print('pure O2', s1, float(s1[2]))
assert s1[0] == 'liquid+gas' and round(float(s1[2]), 2) == 0.41
s2 = state(2); print('air', s2, float(s2[1]))
assert s2[0] == 'liquid' and round(float(s2[1]), -1) == 1870
n = sp.symbols('n', positive=True)
nmin = sp.solve(sp.Eq((Cs + n*CN2)*(Tm - T0)/1000 + Hf + (Cl + n*CN2)*(Tb - Tm)/1000, Q), n)[0]
print('min N2 to avoid vapor', float(nmin)); assert round(float(nmin), 1) == 1.7
print('stages pure O2:', float(Cs*(Tm-T0)/1000), Hf, float(Cl*(Tb-Tm)/1000))
