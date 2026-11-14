# ORG-CBX-051 検証：陰イオン界面活性剤（セッケン・合成洗剤）と非イオン界面活性剤（グルコシド）の区別
import sys
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
P = 'CCCCCCCC/C=C\\CCCCCCCC(=O)[O-].[Na+]'           # オレイン酸ナトリウム
Qs = 'CCCCCCCCCCCCc1ccc(cc1)S(=O)(=O)[O-].[Na+]'     # ドデシルベンゼンスルホン酸ナトリウム
R = 'CCCCCCCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O'  # オクチル β-グルコシド
for name, s in [('P', P), ('Q', Qs), ('R', R)]:
    m = Chem.MolFromSmiles(s)
    cc = m.HasSubstructMatch(Chem.MolFromSmarts('C=C'))
    charged = any(a.GetFormalCharge() != 0 for a in m.GetAtoms())
    hemiacetal = m.HasSubstructMatch(Chem.MolFromSmarts('[CX4;R]([OX2H])[OX2;R]'))
    print(name, CalcMolFormula(m), 'C=C(非芳香族):', cc, 'イオン:', charged, 'ヘミアセタール:', hemiacetal)
mR = Chem.MolFromSmiles(R)
assert CalcMolFormula(mR) == 'C14H28O6'
assert not mR.HasSubstructMatch(Chem.MolFromSmarts('[CX4;R]([OX2H])[OX2;R]'))   # 1位は OH でない → 開環しない
assert Chem.MolFromSmiles(P).HasSubstructMatch(Chem.MolFromSmarts('C=C'))
assert not Chem.MolFromSmiles(Qs).HasSubstructMatch(Chem.MolFromSmarts('C=C'))  # 芳香環は C=C（脂肪族）に当たらない
eqs = ['C17H33COONa + HCl -> C17H33COOH + NaCl',
       '2C17H33COO^- + Mg^2+ -> (C17H33COO)2Mg',
       'C17H33COOH + Br2 -> C17H33Br2COOH',
       'C14H28O6 + H2O -> C6H12O6 + C8H17OH',
       'C17H33COO^- + H2O <=> C17H33COOH + OH^-']
for e in eqs:
    r = check_equation(e); print(e, r); assert r is None
# 判定表（教科書の標準的事実）
table = {
 'ア フェノールフタレイン': {'P': True, 'Q': False, 'R': False},
 'イ 希塩酸': {'P': True, 'Q': False, 'R': False},
 'ウ MgCl2': {'P': True, 'Q': False, 'R': False},
 'エ 臭素水': {'P': True, 'Q': False, 'R': False},
 'オ 銀鏡': {'P': False, 'Q': False, 'R': False},
 'カ オリーブ油': {'P': False, 'Q': False, 'R': False},   # どれも乳化 → 区別できない
 'キ 酸で加熱後フェーリング': {'P': False, 'Q': False, 'R': True},
}
onlyP = [k[0] for k, v in table.items() if v == {'P': True, 'Q': False, 'R': False}]
none_ = [k[0] for k, v in table.items() if not any(v.values())]
onlyR = [k[0] for k, v in table.items() if v == {'P': False, 'Q': False, 'R': True}]
print(onlyP, none_, onlyR)
assert onlyP == ['ア', 'イ', 'ウ', 'エ'] and none_ == ['オ', 'カ'] and onlyR == ['キ']
print('OK')
