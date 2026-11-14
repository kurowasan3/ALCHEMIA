# 検証スクリプト ORG-ARO-094  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/ORG-ARO-094/item.json'
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
import itertools
# C8H10 の芳香族異性体：ベンゼン＋炭素2個の置換
cands = set()
for s in ['c1ccccc1CC']:
    cands.add(Chem.CanonSmiles(s))
for i, j in itertools.combinations(range(6), 2):
    atoms = ['c']*6
    sm = ''
    ring = []
    for k in range(6):
        ring.append('c(C)' if k in (i, j) else 'c')
    sm = ring[0] + '1' + ''.join(ring[1:]) + '1'
    cands.add(Chem.CanonSmiles(sm))
for c in cands: assert CalcMolFormula(Chem.MolFromSmiles(c)) == 'C8H10'
print('C8H10 aromatic', len(cands), cands); assert len(cands) == 4
ox = {'o': 'OC(=O)c1ccccc1C(=O)O', 'm': 'OC(=O)c1cccc(C(=O)O)c1', 'p': 'OC(=O)c1ccc(C(=O)O)cc1', 'Et': 'OC(=O)c1ccccc1'}
for k, v in ox.items(): print(k, CalcMolFormula(Chem.MolFromSmiles(v)))
H, C, O = 1.0, 12, 16
n_tol = 184e3 / (7*C + 8*H); assert abs(n_tol - 2000) < 1e-9
tpa = n_tol/2*0.90*(8*C + 6*H + 4*O)/1e3; print('TPA kg', tpa)
assert sig(tpa, 3) == 149 and abs(0.90/0.25 - 3.6) < 1e-12

print('ALL OK')
