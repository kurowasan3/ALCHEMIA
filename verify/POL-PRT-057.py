# POL-PRT-057 検証：環状ペプチドの加水分解の水・ケルダール法の窒素から組成と分子量を決める
from fractions import Fraction as F
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors as D
Mw = {'Gly': 75, 'Ala': 89, 'Lys': 146}
mass = F(512, 10)                     # mg
water = (F(620, 10) - mass)/18        # mmol（加水分解で増えた質量）
NH3 = F(200, 100)*1 - F(120, 100)     # 0.0500 mol/L H2SO4 20.0 mL = 2.00 mmol H+，NaOH 0.100×12.0 = 1.20 mmol
NH3 = 2*F(5, 100)*20 - F(1, 10)*12
print('水', water, 'mmol  NH3', NH3, 'mmol')
assert water == F(3, 5) and NH3 == F(4, 5)
# g=a（等量），n m = water，(n+l) m = NH3
lm = NH3 - water; gm = (water - lm)/2
print('Lys', lm, 'Gly', gm, 'Ala', gm)
assert lm == gm == F(1, 5)
# 質量の確認
assert gm*(75 + 89) + lm*146 - 18*water == mass
for n in (3, 6, 9):
    M = mass/(water/n); print('n =', n, '分子量', M)
assert mass/(water/6) == 512
# RDKit で cyclo(Gly-Ala-Lys-Gly-Ala-Lys) の分子量・N 数
cyc = 'O=C1CNC(=O)C(C)NC(=O)C(CCCCN)NC(=O)CNC(=O)C(C)NC(=O)C(CCCCN)N1'
m = Chem.MolFromSmiles(cyc); f = D.CalcMolFormula(m); print(f)
mw = 12*m.GetNumAtoms(onlyExplicit=True)  # dummy
from collections import Counter
cnt = Counter(a.GetSymbol() for a in Chem.AddHs(m).GetAtoms())
M = 12*cnt['C'] + 1*cnt['H'] + 14*cnt['N'] + 16*cnt['O']; print('分子量（整数原子量）', M, ' N', cnt['N'])
assert M == 512 and cnt['N'] == 8
# 鎖状なら水は n-1
print('OK')
