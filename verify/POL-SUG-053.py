# 検証スクリプト（東工大 g12）  python3 -I この.py
import json, re
from sympy import symbols, Rational as R, solve, Matrix
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
ITEM = '/home/claude/chemx/tokodai/out/items/POL-SUG-053/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
assert it['problem'].count('\\item') == len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])) == 3
# 構造から OH の数・分子式（RDKit）
smi = {
 'lactose':'OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O',
 'sucrose':'OC[C@H]1O[C@H](O[C@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
 'fructose':'OC[C@@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O'}
patt = Chem.MolFromSmarts('[OX2H]')
nOH = {}
for k,s in smi.items():
    m = Chem.MolFromSmiles(s); nOH[k] = len(m.GetSubstructMatches(patt))
    print(k, rdMolDescriptors.CalcMolFormula(m), 'OH', nOH[k])
assert nOH == {'lactose':8,'sucrose':8,'fructose':5}
# アセチル化の増加 = C2H2O
assert 2*12+2*1+16 == 42
a,b,c = symbols('a b c')
Cu2O = 2*R(635,10)+16; assert Cu2O == 143
eqs = [a+c - R(572,10)/Cu2O, 8*a+8*b+5*c - R(1890,10)/42, 12*a+12*b+6*c - R(2904,10)/44]
sol = solve(eqs,[a,b,c]); print(sol)
assert sol == {a:R(3,10), b:R(2,10), c:R(1,10)}
# 実験2,3だけ：係数行列の階数
M = Matrix([[8,8,5],[12,12,6]]); print('rank', M.rank()); assert M.rank()==2
s2 = solve(eqs[1:], [c, a]); print(s2)   # c と a+b
assert s2[c] == R(1,10) and (s2[a]+b).simplify() == R(1,2)
print('ALL OK')
