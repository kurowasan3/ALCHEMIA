# ORG-HYD-063 verification: Hueckel allyl system, charges, delocalisation energy, HBr addition products
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
x = sp.symbols('x')
A = sp.Matrix([[0,1,0],[1,0,1],[0,1,0]])   # adjacency; H = alpha*I + beta*A
roots = sp.solve(sp.Matrix([[x,1,0],[1,x,1],[0,1,x]]).det(), x); print('x roots', roots)
assert set(roots) == {0, sp.sqrt(2), -sp.sqrt(2)}
al, be = sp.symbols('alpha beta')
vecs = []
for val, mult, vs in A.eigenvects():          # eigenvalue lambda of A -> E = alpha + lambda*beta
    v = vs[0] / sp.sqrt((vs[0].T*vs[0])[0]); vecs.append((val, sp.simplify(v)))
vecs.sort(key=lambda t: -t[0])               # beta<0 -> largest lambda lowest energy
for val, v in vecs: print('E = alpha +', val, '*beta', list(v))
E = [al + val*be for val, _ in vecs]
c2 = [[sp.simplify(c**2) for c in v] for _, v in vecs]
assert c2[0] == [sp.Rational(1,4), sp.Rational(1,2), sp.Rational(1,4)]
assert c2[1] == [sp.Rational(1,2), 0, sp.Rational(1,2)]
def charges(occ):
    q = [sum(occ[j]*c2[j][r] for j in range(3)) for r in range(3)]
    return [1 - qi for qi in q]
cat = charges([2,0,0]); an = charges([2,2,0]); print('cation', cat, 'anion', an)
assert cat == [sp.Rational(1,2), 0, sp.Rational(1,2)] and an == [-sp.Rational(1,2), 0, -sp.Rational(1,2)]
Epi = sp.expand(2*E[0]); deloc = sp.simplify(Epi - 2*(al+be)); print(Epi, deloc)
assert sp.simplify(deloc - 2*(sp.sqrt(2)-1)*be) == 0
val = float(deloc.subs(be, -70)); print('kJ/mol', val, 2*(1.414-1)*70)
assert round(-val) == 58 and round(2*(1.414-1)*70) == 58
for e in ['CH2=CH-CH=CH2 + HBr -> CH3-CHBr-CH=CH2', 'CH2=CH-CH=CH2 + HBr -> CH3-CH=CH-CH2Br']:
    r = check_equation(e); print(e, r); assert r is None
# products are C4H7Br isomers
try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    for smi in ['CC(Br)C=C', 'CC=CCBr']:
        print(smi, CalcMolFormula(Chem.MolFromSmiles(smi))); assert CalcMolFormula(Chem.MolFromSmiles(smi)) == 'C4H7Br'
except ImportError:
    pass
print('OK')
