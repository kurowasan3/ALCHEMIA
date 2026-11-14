# ORG-CBX-009 検証：マレイン酸のメチルエチルエステルの元素分析・加水分解・けん化量
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, nsimplify
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
H, C, O, Na = R(1), R(12), R(16), R(23)
mX, mCO2, mH2O, M = R(79), R(154), R(45), R(158)
nC = mCO2 / 44; nH = 2 * mH2O / 18; mO = mX - nC * C - nH * H; nO = mO / O
n = mX / M
print('C,H,O per molecule:', nC / n, nH / n, nO / n)
assert (nC / n, nH / n, nO / n) == (7, 10, 4)
X = 'CCOC(=O)/C=C\\C(=O)OC'
assert CalcMolFormula(Chem.MolFromSmiles(X)) == 'C7H10O4' and round(MolWt(Chem.MolFromSmiles(X))) == 158
# 原子量 H=1,C=12,O=16 での分子量
assert 7 * 12 + 10 + 4 * 16 == 158
# 加水分解生成物
for s, f in [('OC(=O)/C=C\\C(=O)O', 'C4H4O4'), ('CO', 'CH4O'), ('CCO', 'C2H6O'), ('O=C1OC(=O)C=C1', 'C4H2O3'), ('OC(=O)CCC(=O)O', 'C4H6O4')]:
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == f, s
# けん化に必要な NaOH
mNaOH = R(158, 100) / 158 * 2 * (Na + O + H)
print('NaOH g =', float(mNaOH)); assert mNaOH == R(8, 10)
for eq in ['2C7H10O4 + 15O2 -> 14CO2 + 10H2O', 'CH3OOCCH=CHCOOC2H5 + 2NaOH -> NaOOCCH=CHCOONa + CH3OH + C2H5OH',
           'C4H4O4 -> C4H2O3 + H2O', 'HOOCCH=CHCOOH + H2 -> HOOCCH2CH2COOH']:
    assert check_equation(eq) is None, eq
# P: エタノール（CH3CH(OH)- 構造 → ヨードホルム陽性），Q: メタノール（陰性）
iodo = Chem.MolFromSmarts('[CH3][CX4;H2,H1][OH]')
assert Chem.MolFromSmiles('CCO').HasSubstructMatch(iodo) and not Chem.MolFromSmiles('CO').HasSubstructMatch(iodo)
print('ORG-CBX-009 OK')
