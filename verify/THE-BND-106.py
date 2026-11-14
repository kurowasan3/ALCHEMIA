# 検証スクリプト THE-BND-106  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-BND-106/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|\$', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式の検査', n, '本 OK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from math import sqrt
H, C, O, Zn, NA = 1.0, 12, 16, 65.4, 6.02e23
assert CalcMolFormula(Chem.MolFromSmiles('[O-]C(=O)c1ccc(cc1)C([O-])=O')) == 'C8H4O4-2'
assert CalcMolFormula(Chem.MolFromSmiles('[O-]C(=O)c1ccc(cc1)-c1ccc(cc1)C([O-])=O')) == 'C14H8O4-2'
assert 6 * 0.5 == 3 and 4*2 - 2 + 3*(-2) == 0
M1 = 4*Zn + O + 3*(8*C + 4*H + 4*O); M2 = 4*Zn + O + 3*(14*C + 8*H + 4*O)
r1 = M1/NA/(1.29e-7)**3; r2 = M2/NA/(1.72e-7)**3
print(M1, r1, M2, r2, 0.65/r2, sqrt(3)/2*1.72, sqrt(3)/2*1.29)
assert abs(M1 - 769.6) < 1e-6 and abs(M2 - 997.6) < 1e-6
assert sig(r1, 2) == 0.60 and sig(r2, 2) == 0.33 and round(0.65/r2) == 2
assert sig(sqrt(3)/2*1.72, 3) == 1.49 and sig(sqrt(3)/2*1.29, 3) == 1.12

print('ALL OK')
