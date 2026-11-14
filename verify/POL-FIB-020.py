# POL-FIB-020 検証：縮合重合の分子数（結合 1 つで分子が 1 つ減る），量の不均衡・末端封止・環化
import sys, random; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, simplify, solve, Rational as R, N
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
Mac, Mam, Mw, Mace = 146, 116, 18, 60
assert round(MolWt(Chem.MolFromSmiles('OC(=O)CCCCC(=O)O'))) == 146
assert round(MolWt(Chem.MolFromSmiles('NCCCCCCN'))) == 116
# (1)
nA, nB = R('0.980'), R('1.000')
bonds = 2*nA; mols = nA + nB - bonds; mass = nA*Mac + nB*Mam - bonds*Mw
print('(1) molecules', mols, 'mass', N(mass), 'Mn', N(mass/mols), 'NH2 ends', 2*nB - 2*nA)
assert mols == R('0.02') and round(float(mass/mols), -2) == 11200
# (2) 一般式
r, p, NB = symbols('r p N_B', positive=True)
NA = r*NB
Xn = ((NA + NB)/2) / ((NA + NB)/2 - NA*p)
assert simplify(Xn - (1 + r)/(1 + r - 2*r*p)) == 0
X1 = (nA + nB)/mols; print('Xn(1) =', X1)
pp = solve(1/(1 - p) - X1, p)[0]; print('p needed =', N(pp)); assert round(float(pp), 3) == 0.990
assert simplify(Xn.subs(p, 1) - (1 + r)/(1 - r)) == 0
# (3)
ends = R('2.0e-3') + R('2.0e-3'); total = R('2.5e-3')
lin = ends/2; cyc = total - lin; print('(3) cyclic fraction', N(cyc/total), 'Mn', N(100/total))
assert cyc/total == R(1, 5) and 100/total == 40000
# (5) 酢酸による末端封止
na, nb, nc = R(1), R(1), R('0.0100')
bonds = 2*nb  # NH2 がすべて反応（COOH 2.010 > NH2 2.000）
mols = na + nb + nc - bonds; mass = na*Mac + nb*Mam + nc*Mace - bonds*Mw
print('(5) molecules', mols, 'Mn', N(mass/mols)); assert round(float(mass/mols), -2) == 22700
# 乱数シミュレーション：環をつくらない縮合では分子数＝初期分子数−結合数（union-find）
random.seed(1)
def sim(nacid, namine, frac):
    parent = list(range(nacid + namine))
    def f(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    acid_g = [i for i in range(nacid) for _ in range(2)]
    am_g = [nacid + j for j in range(namine) for _ in range(2)]
    random.shuffle(acid_g); random.shuffle(am_g)
    nb = 0
    target = int(frac*len(acid_g))
    while nb < target and acid_g:
        a = acid_g.pop(); 
        # 環をつくらない相手を探す
        for k in range(len(am_g)):
            if f(am_g[k]) != f(a):
                b = am_g.pop(k); parent[f(a)] = f(b); nb += 1; break
        else:
            break
    comps = len({f(i) for i in range(nacid + namine)})
    return comps, nacid + namine - nb
c1, c2 = sim(490, 500, 1.0); print('sim', c1, c2); assert c1 == c2
for e in ['C6H10O4 + C6H16N2 -> C12H24N2O3 + H2O']:
    assert check_equation(e) is None, e
print('OK')
