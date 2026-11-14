# POL-RES-023 検証：平均官能度（カロザースのモデル）によるフェノール樹脂の網目化の条件
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import Rational as Q, symbols, solve, N, simplify
from chemlib import check_equation
r = symbols('r', positive=True)
def fav(rv):
    rv = Q(rv)
    lim = min(3, 2*rv)          # 不足する側の官能基の数（フェノール1 molあたり）
    return 2*lim/(1+rv)
# 網目化の範囲 1<r<2
for rv in [Q(1,2), Q(4,5), 1, Q(6,5), Q(3,2), Q(9,5), 2, Q(5,2)]:
    print('r=', rv, 'f_av=', fav(rv), '>2?', fav(rv) > 2)
assert fav(Q(4,5)) < 2 and fav(1) == 2 and fav(2) == 2 and fav(Q(3,2)) > 2
print('境界(r<1.5):', solve(4*r/(1+r)-2, r), ' (r>1.5):', solve(6/(1+r)-2, r))
# (2) r=0.8 ノボラック：分子数・平均分子量
n0 = 1 + Q(4,5); bonds = 2*Q(4,5); n_poly = n0 - bonds
mass = 94 + Q(4,5)*30 - Q(4,5)*18
M = mass/n_poly
print('(2) 残る分子', n_poly, ' 質量', mass, ' 平均分子量', M, ' 1分子あたりフェノール', 1/n_poly)
assert M == 518 and 1/n_poly == 5
# 構造で確認：フェノール5 + CH2 橋4
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
pent = 'Oc1ccccc1Cc1cccc(Cc2cccc(Cc3cccc(Cc4ccccc4O)c3O)c2O)c1O'
m = Chem.MolFromSmiles(pent)
print('5量体', CalcMolFormula(m), round(MolWt(m),1))
# 整数原子量で
from collections import Counter
c = Counter(a.GetSymbol() for a in Chem.AddHs(m).GetAtoms())
print(' 整数原子量', c['C']*12 + c['H']*1 + c['O']*16)
assert c['C']*12 + c['H']*1 + c['O']*16 == 518
# (3) ヘキサメチレンテトラミン：不足 0.2 mol CH2 / フェノール1 mol
need_CH2 = 1 - Q(4,5)
hexa = need_CH2/6
print('(3) ヘキサミン mol/フェノール mol', hexa, ' g/(ノボラック100 g)', N(hexa*140*100/mass, 4))
assert check_equation('6HCHO + 4NH3 -> C6H12N4 + 6H2O') is None
# (4) r=1.5, r=1.2 のゲル化点
for rv in [Q(3,2), Q(6,5)]:
    print('(4) r=', rv, 'f_av=', fav(rv), ' p_c=', 2/fav(rv), N(2/fav(rv), 3))
assert 2/fav(Q(3,2)) == Q(5,6)
print('OK')
