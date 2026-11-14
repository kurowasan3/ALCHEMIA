# ORG-ALC-029 検証：グリセリンとベンズアルデヒドの環状アセタールの異性体，メチル化・加水分解生成物の不斉炭素
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def canon(m): return Chem.MolToSmiles(m)
def mirror(smi):
    return canon(Chem.MolFromSmiles(smi.replace('@@', '!').replace('@', '@@').replace('!', '@')))

glycerol = Chem.MolFromSmiles('OCC(O)CO')
# グリセリンの OH 2 個を選んで PhCH< で環にする（五員環・六員環だけ）
ring_smis = {'1,2（五員環）': 'OCC1COC(c2ccccc2)O1', '1,3（六員環）': 'OC1COC(c2ccccc2)OC1'}
total = 0; chiral_total = 0
for k, s in ring_smis.items():
    m = Chem.MolFromSmiles(s)
    assert rdMolDescriptors.CalcMolFormula(m) == 'C10H12O3'
    opts = StereoEnumerationOptions(unique=True)
    isos = sorted({canon(x) for x in EnumerateStereoisomers(m, options=opts)})
    chiral = [x for x in isos if mirror(x) != x]
    print(k, '立体異性体', len(isos), isos, ' キラル', len(chiral))
    total += len(isos); chiral_total += len(chiral)
    # 環の大きさ
    ri = m.GetRingInfo(); print('  環の員数', sorted(len(r) for r in ri.AtomRings()))
assert total == 6 and chiral_total == 4
# 六員環の cis/trans は互いに鏡像でない（ジアステレオマー）・それぞれアキラル
six = sorted({canon(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(ring_smis['1,3（六員環）']))})
assert len(six) == 2 and all(mirror(x) == x for x in six)
# 分子量
print('M(C10H12O3)=', round(Descriptors.MolWt(Chem.MolFromSmiles('OCC1COC(c2ccccc2)O1')), 2))
# メチル化（遊離 OH → OCH3）→ 酸加水分解（アセタール → PhCHO + ジオール）
p5 = Chem.MolFromSmiles('COCC(O)CO')   # 五員環から：3-O-メチルグリセリン
p6 = Chem.MolFromSmiles('OCC(OC)CO')   # 六員環から：2-O-メチルグリセリン
for nm, m in (('五員環→', p5), ('六員環→', p6)):
    cc = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    print(nm, Chem.MolToSmiles(m), rdMolDescriptors.CalcMolFormula(m), '不斉炭素', len(cc))
assert len(Chem.FindMolChiralCenters(p5, includeUnassigned=True, useLegacyImplementation=False)) == 1
assert len(Chem.FindMolChiralCenters(p6, includeUnassigned=True, useLegacyImplementation=False)) == 0
assert rdMolDescriptors.CalcMolFormula(p5) == rdMolDescriptors.CalcMolFormula(p6) == 'C4H10O3'
# 反応式
for eq in ['C3H8O3 + C7H6O -> C10H12O3 + H2O', 'C11H14O3 + H2O -> C7H6O + C4H10O3']:
    assert check_equation(eq) is None, eq
# 質量計算：アセタール 9.0 g（0.050 mol）から PhCHO
n = 9.0 / 180.0; print('PhCHO g', n * 106)
print('ALL OK')
