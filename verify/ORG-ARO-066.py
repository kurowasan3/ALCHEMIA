# ORG-ARO-066 検証：2-ブタノンとフェノール2分子の縮合（ビスフェノール）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import FindMolChiralCenters
import sympy as sp
# X の候補：C4H8O のカルボニル化合物（ブタナール，2-メチルプロパナール，2-ブタノン）
cands = {'ブタナール': 'CCCC=O', '2-メチルプロパナール': 'CC(C)C=O', '2-ブタノン': 'CCC(C)=O'}
def iodoform(smi):  # CH3-CO- をもつか
    return Chem.MolFromSmiles(smi).HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)'))
def aldehyde(smi):
    return Chem.MolFromSmiles(smi).HasSubstructMatch(Chem.MolFromSmarts('[CX3H1](=O)'))
ok = [n for n, s in cands.items() if iodoform(s) and not aldehyde(s)]
print('X =', ok); assert ok == ['2-ブタノン']
Y = Chem.MolFromSmiles('CCC(C)(c1ccc(O)cc1)c1ccc(O)cc1')
print('Y formula', CalcMolFormula(Y)); assert CalcMolFormula(Y) == 'C16H18O2'
cc = FindMolChiralCenters(Y, includeUnassigned=True, useLegacyImplementation=False)
print('chiral centers', cc); assert cc == []
assert check_equation('C4H8O + 2C6H6O -> C16H18O2 + H2O') is None
H, C, O = 1.0, 12, 16
MX = 4*C + 8*H + O; MY = 16*C + 18*H + 2*O
print('MX, MY', MX, MY)
mY = sp.Rational(72, 10) / sp.nsimplify(MX) * sp.Rational(8, 10) * sp.nsimplify(MY)
print('Y の質量', float(mY), 'g -> 19 g'); assert round(float(mY)) == 19
# フェノール性 OH 2 個 → NaOH 2 mol
print('phenolic OH', len(Y.GetSubstructMatches(Chem.MolFromSmarts('c[OX2H]'))))
