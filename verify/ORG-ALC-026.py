# ORG-ALC-026 検証：アセトアルデヒドとエチレングリコールの環状アセタール P，異性体との区別，銀の質量
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
P = Chem.MolFromSmiles('CC1OCCO1')            # 2-メチル-1,3-ジオキソラン
EA = Chem.MolFromSmiles('CCOC(C)=O')          # 酢酸エチル
for m in (P, EA): assert rdMolDescriptors.CalcMolFormula(m) == 'C4H8O2'
print('M(P)=', round(Descriptors.MolWt(P), 1))
# P はアセタール炭素（O 2 個と結合，OH なし），アルデヒド基なし，エステルなし
acetal = Chem.MolFromSmarts('[CX4]([OX2][#6])[OX2][#6]')
assert P.HasSubstructMatch(acetal) and not P.HasSubstructMatch(Chem.MolFromSmarts('[CX3H1]=O'))
assert not P.HasSubstructMatch(Chem.MolFromSmarts('C(=O)O'))
ring = P.GetRingInfo().AtomRings(); assert [len(r) for r in ring] == [5]
# 反応式
eqs = ['CH3CHO + HOCH2CH2OH -> C4H8O2 + H2O',
       'C4H8O2 + H2O -> CH3CHO + HOCH2CH2OH',
       'CH3COOC2H5 + NaOH -> CH3COONa + C2H5OH',
       'CH3CHO + 2[Ag(NH3)2]+ + 3OH- -> CH3COO- + 2Ag + 4NH3 + 2H2O']
for e in eqs:
    r = check_equation(e); print(e, r); assert r is None
# 8.8 g の P → CH3CHO 0.10 mol → Ag 0.20 mol × 108
n = 8.8 / 88; mAg = 2 * n * 108; print('Ag g', mAg); assert round(mAg) == 22
print('ALL OK')
