# 検証スクリプト POL-RES-031  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/POL-RES-031/item.json'
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

from sympy import Rational, symbols, solve, nsimplify
from rdkit import Chem
H, C = 1.0, 12.0
ne = 560 / (2*C + 4*H); nh = 21.0 / (6*C + 12*H); nb = 21.0 / (4*C + 8*H)
assert abs(ne - 20.0) < 1e-9 and abs(nh - 0.25) < 1e-9 and abs(nb - 0.375) < 1e-9
b2 = nh / (2*(ne + nh)) * 1000; b3 = nb / (2*(ne + nb)) * 1000
print('(2)', b2, '(3)', b3)
assert sig(b2, 2) == 6.2 and sig(b3, 2) == 9.2 and b3 > b2
x = symbols('x', positive=True)
sol = solve(x / (2*(20 + x)) * 1000 - 10, x)[0]
m = float(sol) * 84
print('(4)', float(sol), m)
assert sig(m, 3) == 34.3
# 部分構造：主鎖 C-C(ブチル)-C-C
mol = Chem.MolFromSmiles('CCC(CCCC)CC')
atom = [a for a in mol.GetAtoms() if a.GetDegree() == 3][0]
assert atom.GetTotalNumHs() == 1
assert '34.3' in it['answer'] and '6.2' in it['answer'] and '9.2' in it['answer']

print('ALL OK')
