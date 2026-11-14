# ORG-HYD-048 検証：アルケン混合物のオゾン分解生成物をヨードホルム反応と銀鏡反応で定量
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as R, N
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def ozonolysis(smi):
    m = Chem.MolFromSmiles(smi)
    rxn_out = []
    from rdkit.Chem import AllChem
    rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]=O.[C:2]=O')
    ps = rxn.RunReactants((m,))
    prods = set()
    for p in ps:
        s = []
        for q in p:
            Chem.SanitizeMol(q); s.append(Chem.MolToSmiles(q))
        prods.add(tuple(sorted(s)))
    assert len(prods) == 1
    return list(prods)[0]
iodo = Chem.MolFromSmarts('[CH3]C(=O)[#1,#6]')
iodo2 = Chem.MolFromSmarts('[CH3][CX3](=O)')
ald = Chem.MolFromSmarts('[CX3H1](=O)')
def counts(smis):
    ci = sum(len({m[1] for m in Chem.MolFromSmiles(s).GetSubstructMatches(iodo2)}) for s in smis)  # CH3-C(=O) の C=O 炭素 1 個につき CHI3 1 個
    ag = 2*sum(len(Chem.MolFromSmiles(s).GetSubstructMatches(ald)) for s in smis)
    return ci, ag
alk = {'2-methyl-2-pentene': 'CC(C)=CCC', '3-methyl-2-pentene': 'CC=C(C)CC', '2,3-dimethyl-2-butene': 'CC(C)=C(C)C'}
res = {}
for k, s in alk.items():
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C6H12'
    pr = ozonolysis(s); res[k] = counts(pr); print(k, pr, 'CHI3, Ag per mol =', res[k])
assert res['2-methyl-2-pentene'] == (1, 2) and res['3-methyl-2-pentene'] == (2, 2) and res['2,3-dimethyl-2-butene'] == (2, 0)
x, y = symbols('x y')
MCHI3 = 12 + 1 + 3*127; assert MCHI3 == 394
nI = R('23.6')/MCHI3; nAg = R('10.8')/108
print('n CHI3 =', N(nI), 'n Ag =', N(nAg))
sol = solve([x + y - R('0.0500'), x + 2*y - R('0.0600')], [x, y])
print(sol); assert sol[x] == R('0.04') and sol[y] == R('0.01')
assert abs(float(nI) - 0.0600) < 0.0002 and nAg == R('0.1') and 2*(sol[x] + sol[y]) == nAg
print('mole % of 3-methyl-2-pentene =', N(sol[y]/R('0.05')*100))
for e in ['CH3COCH3 + 3I2 + 4NaOH -> CHI3 + CH3COONa + 3NaI + 3H2O',
          'CH3CHO + 3I2 + 4NaOH -> CHI3 + HCOONa + 3NaI + 3H2O',
          'C2H5CHO + 2[Ag(NH3)2]+ + 3OH- -> C2H5COO- + 2Ag + 4NH3 + 2H2O']:
    assert check_equation(e) is None, (e, check_equation(e))
print('OK')
