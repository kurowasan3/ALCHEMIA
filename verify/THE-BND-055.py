# THE-BND-055 検証：価電子の総数の偶奇と不対電子・等電子
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
val = {'C': 4, 'N': 5, 'O': 6, 'Cl': 7}
sp = {'CO': (['C','O'], 0), 'NO': (['N','O'], 0), 'NO2': (['N','O','O'], 0), 'CN-': (['C','N'], -1),
      'ClO2': (['Cl','O','O'], 0), 'N2O': (['N','N','O'], 0), 'NO+': (['N','O'], +1), 'N2': (['N','N'], 0)}
cnt = {k: sum(val[a] for a in at) - q for k, (at, q) in sp.items()}
print(cnt)
assert cnt == {'CO': 10, 'NO': 11, 'NO2': 17, 'CN-': 10, 'ClO2': 19, 'N2O': 16, 'NO+': 10, 'N2': 10}
odd = [k for k, v in cnt.items() if v % 2]
print('奇数:', odd); assert odd == ['NO', 'NO2', 'ClO2']
iso = [k for k, v in cnt.items() if v == 10 and len(sp[k][0]) == 2 and k != 'N2']
print('N2 と等電子（2原子・10個）:', iso); assert iso == ['CO', 'CN-', 'NO+']
# CO：三重結合（6電子）＋非共有電子対 2 組（4電子）＝10
assert 3*2 + 2*2 == cnt['CO']
# RDKit でラジカル電子数を確認
from rdkit import Chem
for smi, name in [('[N]=O', 'NO'), ('[O]N=O', 'NO2'), ('[C-]#[O+]', 'CO'), ('[C-]#N', 'CN-'), ('N#[O+]', 'NO+')]:
    m = Chem.MolFromSmiles(smi)
    print(name, 'ラジカル電子', sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()))
for e in ['2NO2 <=> N2O4', '2NO + O2 -> 2NO2']:
    assert check_equation(e) is None, e
print('OK')
