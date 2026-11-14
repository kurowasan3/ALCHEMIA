# POL-SUG-033 検証：ペンタアセチルグルコースの収率，1位だけを H にした化合物 X と鎖状の還元生成物 Y
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import sympy as sp
H, C, O = 1, 12, 16
M = lambda c, h, o: c * C + h * H + o * O
Mglc, Mpa, MX, MY = M(6, 12, 6), M(16, 22, 11), M(6, 12, 5), M(6, 14, 6)
print('M', Mglc, Mpa, MX, MY); assert (Mglc, Mpa, MX, MY) == (180, 390, 164, 182)
n = sp.Rational(9, 1) / Mglc                      # 9.0 g
theo = n * Mpa; yld = sp.Rational(156, 10) / theo
print('理論量', theo, 'g  収率', float(yld)); assert theo == sp.Rational(39, 2) and yld == sp.Rational(4, 5)
nX = sp.Rational(156, 10) / Mpa; mX = nX * MX
print('X の質量', float(mX)); assert round(float(mX), 1) == 6.6
# 構造の確認（立体は問わない）
glc = Chem.MolFromSmiles('OCC1OC(O)C(O)C(O)C1O')
pa = Chem.MolFromSmiles('CC(=O)OCC1OC(OC(C)=O)C(OC(C)=O)C(OC(C)=O)C1OC(C)=O')
X = Chem.MolFromSmiles('OCC1OCC(O)C(O)C1O')       # 1,5-アンヒドログルシトール
Y = Chem.MolFromSmiles('OCC(O)C(O)C(O)C(O)CO')    # グルシトール（ソルビトール）
for nm, m, f in (('glc', glc, 'C6H12O6'), ('pa', pa, 'C16H22O11'), ('X', X, 'C6H12O5'), ('Y', Y, 'C6H14O6')):
    assert rdMolDescriptors.CalcMolFormula(m) == f, nm
# X にはヘミアセタール（C に OH と環の O が同時に結合）もアルデヒドもない
hemi = Chem.MolFromSmarts('[CX4]([OX2H])[OX2][#6]')
cho = Chem.MolFromSmarts('[CX3H1]=O')
assert glc.HasSubstructMatch(hemi) and not X.HasSubstructMatch(hemi) and not X.HasSubstructMatch(cho)
assert not Y.HasSubstructMatch(hemi) and not Y.HasSubstructMatch(cho)
# 酢酸エステルの数 = 5
assert len(pa.GetSubstructMatches(Chem.MolFromSmarts('OC(C)=O'))) == 5
for eq in ['C6H12O6 + 5(CH3CO)2O -> C16H22O11 + 5CH3COOH',
           'C16H22O11 + HBr -> C14H19BrO9 + CH3COOH',
           'C14H19BrO9 + H2 -> C14H20O9 + HBr',
           'C14H20O9 + 4H2O -> C6H12O5 + 4CH3COOH']:
    r = check_equation(eq); print(eq, r); assert r is None
print('ALL OK')
