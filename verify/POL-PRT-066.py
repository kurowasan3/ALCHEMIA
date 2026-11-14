# POL-PRT-066 元素分析と加水分解生成物によるアミド・エステル化合物の構造決定
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import FindMolChiralCenters
for e in ['CH3CONHCH2COOC2H5 + 2H2O -> CH3COOH + H2NCH2COOH + C2H5OH',
          'CH3CONHCH2COOC2H5 + 2NaOH -> CH3COONa + H2NCH2COONa + C2H5OH',
          '4C6H11NO3 + 29O2 -> 24CO2 + 22H2O + 2N2']:
    r = check_equation(e); assert r is None, (e, r)
# 元素分析（H=1, C=12, N=14, O=16）
w = R(290, 10); mC = R(528, 10)*12/44; mH = R(198, 10)*2/18; mN = R(280, 100)
mO = w - mC - mH - mN
nC, nH, nN, nO = mC/12, mH/1, mN/14, mO/16
print('mmol C H N O', nC, nH, nN, nO)
ratio = [x/nN for x in (nC, nH, nN, nO)]; print('ratio', ratio)
assert ratio == [6, 11, 1, 3]
Mf = 6*12 + 11 + 14 + 48; print('M =', Mf); assert Mf == 145 and 2*Mf > 200
A = Chem.MolFromSmiles('CC(=O)NCC(=O)OCC')
assert CalcMolFormula(A) == 'C6H11NO3'
assert len(FindMolChiralCenters(A, includeUnassigned=True, useLegacyImplementation=False)) == 0
gly = Chem.MolFromSmiles('NCC(=O)O'); assert len(FindMolChiralCenters(gly, includeUnassigned=True, useLegacyImplementation=False)) == 0
# 生成物の和 − 2H2O = A
print('A formula OK; all OK')
