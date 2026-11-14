# POL-SUG-040 検証：マルトース・スクロース混合物と Cu2O の量
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
H, C, O, Cu = 1.0, 12.0, 16.0, 63.5
M12 = 12*C + 22*H + 11*O; Mcu = 2*Cu + O
print('C12H22O11', M12, 'Cu2O', Mcu)
n = 6.84/M12; x = 4*0.2145/Mcu; y = n - x   # 操作1 は S の 1/4
print('total', n, 'maltose', x, 'sucrose', y, 'mass% maltose', x/n*100)
op2 = (x+2*y)/2*Mcu; op3 = 2*n/4*Mcu   # 操作2 は 1/2，操作3 は 1/4
print('op2', op2, 'op3', op3)
assert round(x/n*100, 1) == 30.0
assert round(op2, 2) == 2.43 and round(op3, 2) == 1.43
# 完全加水分解後だけでは組成によらない
for xm in (0.0, 0.006, 0.02):
    assert abs(2*((xm)+(n-xm))*Mcu - 2*n*Mcu) < 1e-12
sucrose = 'OC[C@H]1O[C@H](O[C@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O'
maltose = 'OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O'
for s in (sucrose, maltose):
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C12H22O11'
for e in ['C12H22O11 + H2O -> 2C6H12O6']:
    assert check_equation(e) is None
print('OK')
