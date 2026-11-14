# 検証スクリプト POL-RUB-011  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/POL-RUB-011/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, formula_mass
from sympy import symbols, solve, Rational as R
from rdkit import Chem
it = json.load(open(ITEM, encoding='utf-8'))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 4, (n_items, n_ans)
eqs = sorted({b for k in ('problem', 'answer', 'solution') for b in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[k])
              if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm|\$', b)})
print(eqs); assert eqs == ['HO(CH2)4OH -> CH2=CH-CH=CH2 + 2H2O']
for e in eqs: assert check_equation(e) is None, check_equation(e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
M = {'H': 1.0, 'C': 12, 'N': 14, 'Br': 80}
an, bd = formula_mass('C3H3N', M), formula_mass('C4H6', M); assert (an, bd) == (53, 54)
x = symbols('x'); xs = solve(R(14) / (53 + 54*x) * 100 - R(651, 100), x)[0]; print('x', float(xs)); assert round(float(xs), 2) == 3.00
assert sig(14/215*100, 3) == 6.51
br = 10.0/215*3*formula_mass('Br2', M); print('Br2', br); assert sig(br, 3) == 22.3
# 1,4-単位と 1,2-単位のモデル（主鎖の切れ目を H にしたもの）の C=C 数
for smi in ('CC=CC', 'CC(C=C)'):
    m = Chem.MolFromSmiles(smi)
    n = sum(1 for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE); print(smi, n); assert n == 1
print('全検証OK:', ITEM)
