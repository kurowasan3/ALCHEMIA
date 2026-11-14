# POL-SUG-047 検証：結合位置（1,4／1,6／1,3）によるグルコース単位の過ヨウ素酸酸化
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
from fractions import Fraction as F
# 鎖中のグルコース単位をメチル（隣の単位の代わり）で表したモデル：O1 と結合位置の O をメチル化
models = {'1,4': 'COC1OC(CO)C(OC)C(O)C1O', '1,6': 'COC1OC(COC)C(O)C(O)C1O', '1,3': 'COC1OC(CO)C(O)C(OC)C1O'}
def vicinal(m):   # OH をもつ炭素どうしの C-C 結合（隣接ジオール）
    oh = {a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'C' and any(n.GetSymbol() == 'O' and n.GetTotalNumHs() == 1 for n in a.GetNeighbors())}
    return [b for b in m.GetBonds() if b.GetBeginAtomIdx() in oh and b.GetEndAtomIdx() in oh]
cnt = {}
for k, s in models.items():
    m = Chem.MolFromSmiles(s); v = vicinal(m)
    # 中央の炭素（両側が切れる炭素）の数＝ギ酸
    atoms = [b.GetBeginAtomIdx() for b in v] + [b.GetEndAtomIdx() for b in v]
    hcooh = sum(1 for a in set(atoms) if atoms.count(a) == 2)
    cnt[k] = (len(v), hcooh); print(k, '切れる結合', len(v), 'ギ酸', hcooh)
assert cnt == {'1,4': (1, 0), '1,6': (2, 1), '1,3': (0, 0)}
# 単位の式量：C6H10O5=162 → 1,4：+O−H2O → C6H8O5=160；1,6：+2O−H2O−HCOOH → C5H6O4=130
M = lambda C, H, O: 12*C + H + 16*O
assert M(6, 10, 5) == 162 and M(6, 8, 5) == 160 and M(5, 6, 4) == 130
n = F(162, 10)/162   # 16.2 g → 0.100 mol
print('アミロース 生成物', n*160, 'g  IO4-', n*1, '  デキストラン 生成物', n*130, 'g  IO4-', n*2, ' HCOOH', n)
# 混合：IO4- 0.130 mol → 1,6 の割合
x = (F(130, 1000) - n)/n; print('1,6 の割合', x, ' HCOOH', x*n)
assert x == F(3, 10)
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['C6H10O5 + IO4- -> C6H8O5 + IO3- + H2O', 'C6H10O5 + 2IO4- -> C5H6O4 + HCOOH + 2IO3- + H2O']:
    r = check_equation(e); print(e, r); assert r is None
print('OK')
