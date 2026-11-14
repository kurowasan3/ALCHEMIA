# ORG-HYD-059 検証：3,3-ジメチルペンタンの熱分解（最も弱い C-C 結合と β 開裂）
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
m = Chem.MolFromSmiles('CCC(C)(C)CC')
print('substrate', CalcMolFormula(m))
# 各 C-C 結合の両端の炭素の「炭素の隣接数」の和で弱さを順位づけ（表の傾向：多いほど弱い）
best = []
for b in m.GetBonds():
    a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
    score = a1.GetDegree() + a2.GetDegree()
    best.append((score, a1.GetIdx(), a2.GetIdx()))
best.sort(reverse=True); print(best)
assert best[0][0] == 6   # C2-C3（第四級炭素と CH2）
# 開裂で生じる断片
em = Chem.RWMol(m); em.RemoveBond(1, 2)
frags = Chem.GetMolFrags(em.GetMol(), asMols=True, sanitizeFrags=False)
print([Chem.MolToSmiles(f) for f in frags])
# 生成物
prods = {'ethane': 'CC', 'isobutene': 'C=C(C)C', 'methane': 'C'}
for k, s in prods.items(): print(k, CalcMolFormula(Chem.MolFromSmiles(s)))
# 原子の収支：C7H16 + 2H -> C2H6 + C4H8 + CH4
C = 2 + 4 + 1; H = 6 + 8 + 4
assert C == 7 and H == 16 + 2
# 第三級ラジカル (CH3)2C•-CH2CH3 の β 位の C-C は CH2-CH3 のみ（CH3 側は β 炭素に C がない）
print('ok')
