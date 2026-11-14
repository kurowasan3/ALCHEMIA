# ORG-DET-035 verification: identification truth table and equations
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from itertools import combinations, permutations
eqs = ['Na2SO3 + 2HCl -> 2NaCl + SO2 + H2O', 'SO3^2- + I2 + H2O -> SO4^2- + 2I- + 2H+',
       'C6H4(OH)COOH + NaHCO3 -> C6H4(OH)COONa + H2O + CO2', 'Na2CO3 + HCl -> NaCl + NaHCO3',
       'NaHCO3 + HCl -> NaCl + H2O + CO2', 'C6H5OH + Na2CO3 -> C6H5ONa + NaHCO3',
       'C2H5OH + 4I2 + 6NaOH -> CHI3 + HCOONa + 5NaI + 5H2O']
for e in eqs:
    r = check_equation(e); print(e, r); assert r is None
A = ['phenol', 'salicylic', 'acetic', 'ethanol']; B = ['Na2SO4', 'Na2SO3', 'Na2CO3', 'NaHCO3']
pKa = {'salicylic': 3.0, 'acetic': 4.76, 'phenol': 10.0, 'ethanol': 16}
H2CO3_pKa1 = 6.4
def hcl(x):  # (gas immediately on first drop, odorless?, pungent?)
    return {'NaHCO3': ('immediate', 'odorless'), 'Na2CO3': ('delayed', 'odorless'),
            'Na2SO3': ('gas', 'pungent')}.get(x, (None, None))
def iodoform(x): return x == 'ethanol'
def fecl3_violet(x): return x in ('phenol', 'salicylic')
def gas_with(y, x):  # y added to x
    if x == 'NaHCO3': return pKa.get(y, 99) < H2CO3_pKa1
    return False
sols = []
for a in combinations(A, 2):
    for b in combinations(B, 2):
        for W, X, Y, Z in permutations(a + b):
            ok = (hcl(X) == ('immediate', 'odorless') and hcl(Z)[1] == 'pungent'
                  and hcl(W) == (None, None) and hcl(Y) == (None, None)
                  and Z == 'Na2SO3' and iodoform(W) and fecl3_violet(Y) and gas_with(Y, X))
            if ok: sols.append((W, X, Y, Z))
print(sols); assert sols == [('ethanol', 'NaHCO3', 'salicylic', 'Na2SO3')]
print('OK')
