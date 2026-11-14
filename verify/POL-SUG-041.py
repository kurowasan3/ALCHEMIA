# POL-SUG-041 検証：アルドースの環の大きさ・不斉炭素・還元性
import sys
from rdkit import Chem
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def nchiral(smi):
    m = Chem.MolFromSmiles(smi)
    return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))

# 環の員数 = 付加する OH の炭素番号 k に対し (C1..Ck) + O = k+1
def ring_sizes(ncarbon, oh_positions):
    return sorted(k + 1 for k in oh_positions if k >= 2 and k + 1 in (5, 6))
cases = {
    'グリセルアルデヒド': (3, [2, 3]),
    'エリトロース': (4, [2, 3, 4]),
    'リボース': (5, [2, 3, 4, 5]),
    '2-デオキシリボース': (5, [3, 4, 5]),
    'グルコース': (6, [2, 3, 4, 5, 6]),
}
res = {k: ring_sizes(*v) for k, v in cases.items()}
print(res)
assert res == {'グリセルアルデヒド': [], 'エリトロース': [5], 'リボース': [5, 6],
               '2-デオキシリボース': [5, 6], 'グルコース': [5, 6]}

# 環状構造をRDKitで作り，環の員数と分子式を確認
smi = {
    'glc_furanose': 'OCC(O)C1OC(O)C(O)C1O',
    'glc_pyranose': 'OCC1OC(O)C(O)C(O)C1O',
    'glc_chain': 'OCC(O)C(O)C(O)C(O)C=O',
    'rib_furanose': 'OCC1OC(O)C(O)C1O',
    'rib_pyranose': 'OC1COC(O)C(O)C1O',
    'rib_chain': 'OCC(O)C(O)C(O)C=O',
    'ery_furanose': 'OC1COC(O)C1O',
    'dri_furanose': 'OCC1OC(O)CC1O',
    'dri_pyranose': 'OC1COC(O)CC1O',
}
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
for k, s in smi.items():
    m = Chem.MolFromSmiles(s)
    rings = [len(r) for r in m.GetRingInfo().AtomRings()]
    print(k, CalcMolFormula(m), rings, 'chiral C =', nchiral(s))
assert CalcMolFormula(Chem.MolFromSmiles(smi['glc_furanose'])) == 'C6H12O6'
assert nchiral(smi['glc_furanose']) == 5 and nchiral(smi['glc_pyranose']) == 5
assert nchiral(smi['glc_chain']) == 4
# 五員環グルコースの環外側鎖：C4 に結合した CH(OH)CH2OH
# 水中のリボースの形：鎖状 + 五員環(α,β) + 六員環(α,β) = 5
assert 1 + 2 * len(res['リボース']) == 5
# 還元性：1位が OH のまま（ヘミアセタール）または鎖状アルデヒド → 還元性あり
# メチルリボシド（1位 OCH3）はアセタールで開環しない → 還元性なし
print(check_equation('C5H10O5 + CH4O -> C6H12O5 + H2O'))
assert check_equation('C5H10O5 + CH4O -> C6H12O5 + H2O') is None

# 2-デオキシリボース六員環：環の炭素のうち OH をもたないもの（C2 と C5）が 2 個
m = Chem.MolFromSmiles(smi['dri_pyranose'])
ring = m.GetRingInfo().AtomRings()[0]
noOH = []
for idx in ring:
    a = m.GetAtomWithIdx(idx)
    if a.GetSymbol() != 'C':
        continue
    hasOH = any(n.GetSymbol() == 'O' and n.GetTotalNumHs() == 1 for n in a.GetNeighbors())
    nO = sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'O')
    if not hasOH:
        noOH.append((idx, nO))
print('deoxyribopyranose ring C without OH:', noOH)
assert len(noOH) == 2 and all(nO <= 1 for _, nO in noOH)
print('OK')
