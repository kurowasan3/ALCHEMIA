# THE-BND-052 検証：混成軌道・非共有電子対・配位結合による混成の変化（RDKit）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
r = check_equation('BF3 + NH3 -> BF3NH3'); print('BF3+NH3', r or 'OK'); assert r is None
r = check_equation('NH3 + H+ -> NH4+'); assert r is None
cases = {'H2O': ('O', 0, 'SP3'), 'BF3': ('FB(F)F', 1, 'SP2'), 'CO2': ('O=C=O', 1, 'SP'),
         'BF4-': ('F[B-](F)(F)F', 1, 'SP3'), 'NH4+': ('[NH4+]', 0, 'SP3'),
         'F3B-NH3 (B)': ('F[B-](F)(F)[NH3+]', 1, 'SP3'), 'F3B-NH3 (N)': ('F[B-](F)(F)[NH3+]', 4, 'SP3'),
         'NH3': ('N', 0, 'SP3')}
for name, (smi, idx, hyb) in cases.items():
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); a = m.GetAtomWithIdx(idx)
    # 非共有電子対の数 = (価電子 - 形式電荷 - 結合数)/2
    val = {'O': 6, 'B': 3, 'C': 4, 'N': 5}[a.GetSymbol()]
    lp = (val - a.GetFormalCharge() - a.GetTotalValence()) // 2
    print(name, a.GetSymbol(), a.GetHybridization(), 'lone pairs', lp)
    assert str(a.GetHybridization()) == hyb
# 結合角（実測・理想）：CO2 180, BF3 120, NH4+ 109.5, NH3 106.7, H2O 104.5
angles = {'CO2': 180, 'BF3': 120, 'NH4+': 109.5, 'NH3': 106.7, 'H2O': 104.5}
order = sorted(angles, key=angles.get, reverse=True); print(order)
assert order == ['CO2', 'BF3', 'NH4+', 'NH3', 'H2O']
print('OK')
