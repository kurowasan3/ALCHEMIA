# ORG-CBX-055 検証：不飽和カルボン酸の酸化開裂とラクトンの環
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
f = lambda s: CalcMolFormula(Chem.MolFromSmiles(s))
X = 'CC=CCCC(=O)O'; assert f(X) == 'C6H10O2'
# 酸化開裂：C4=C5 を切る
assert f('CC(=O)O') == 'C2H4O2' and f('OC(=O)CCC(=O)O') == 'C4H6O4'
Y1, Y2 = 'CCC(O)CCC(=O)O', 'CC(O)CCCC(=O)O'
assert f(Y1) == f(Y2) == 'C6H12O3'
L1, L2 = 'CCC1CCC(=O)O1', 'CC1CCCC(=O)O1'
for L, n in ((L1, 5), (L2, 6)):
    m = Chem.MolFromSmiles(L); assert f(L) == 'C6H10O2'
    assert [len(r) for r in Chem.GetSymmSSSR(m)] == [n]
# (4) 各ヒドロキシ酸の分子内エステル
for s, n in (('CC1OC1=O', 3), ('CC1CC(=O)O1', 4), ('O=C1CCCO1', 5), ('O=C1CCCCO1', 6)):
    assert [len(r) for r in Chem.GetSymmSSSR(Chem.MolFromSmiles(s))] == [n]
# 他の位置の二重結合では別の開裂生成物（例：CCC=CCC(=O)O → プロピオン酸＋マロン酸）
print('OK')
