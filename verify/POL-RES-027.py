# POL-RES-027 ヒドロキシ酸の縮合で生じる環状物と鎖状物（脱水量の数え方と平均重合度）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, symbols, solve
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert check_equation('2HOCH2COOH -> C4H4O4 + 2H2O') is None
# 原子量 H=1, C=12, O=16 による式量
Mg = 2*12 + 4*1 + 3*16; assert Mg == 76
unit = Mg - 18; assert unit == 58
gly = Chem.MolFromSmiles('O=C1COC(=O)CO1'); print('glycolide', CalcMolFormula(gly))
assert CalcMolFormula(gly) == 'C4H4O4'
# データ
total_A = R(7600, 100); cyc = R(580, 100); chain = R(5328, 100)
water = total_A - cyc - chain
nw = water/18
print('water g, mol', water, nw); assert nw == R(94, 100)
n_cyc = cyc/unit; print('A in rings', n_cyc); assert n_cyc == R(1, 10)
n_chain_units = 1 - n_cyc
n_bonds_chain = nw - n_cyc          # 環: 単位数 = エステル結合数 = 脱水数
N_chain = n_chain_units - n_bonds_chain
print('chains mol', N_chain, 'avg n', n_chain_units/N_chain)
assert N_chain == R(6, 100) and n_chain_units/N_chain == 15
# 鎖状物の質量の整合
assert n_chain_units*Mg - 18*n_bonds_chain == chain
# 環状物がすべて二量体（グリコリド）のとき
print('glycolide mol', n_cyc/2)
print('all OK')
