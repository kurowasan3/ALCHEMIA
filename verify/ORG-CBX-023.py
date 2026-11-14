# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/out/items/ORG-CBX-023/item.json'
it = json.load(open(ITEM, encoding='utf-8'))


def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0:
            return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i + 4:j - 1]); i = j


eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
# 小問数と答えの番号の一致
nq = it['problem'].count('\\item')
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問', nq, '答え', na)


def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import sympy as sp
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
H, C, O = 1.0, 12, 16
M = 4*C + 4*H + 4*O
assert M == 116
n = sp.Rational(116, 10) / 116
assert sig(float(n*22.4), 3) == 2.24
assert sig(float(2*n/1.00*1000), 3) == 200
nw = 0.630 / 18.0
mm = nw * M
pct = mm / 11.6 * 100
print('水', nw, 'マレイン酸', mm, 'g', pct, '%')
assert sig(pct, 3) == 35.0
mal = Chem.MolFromSmiles('OC(=O)/C=C\\C(=O)O'); fum = Chem.MolFromSmiles('OC(=O)/C=C/C(=O)O')
anh = Chem.MolFromSmiles('O=C1OC(=O)C=C1')
assert CalcMolFormula(mal) == CalcMolFormula(fum) == 'C4H4O4'
assert CalcMolFormula(anh) == 'C4H2O3'
assert Chem.MolToSmiles(mal) != Chem.MolToSmiles(fum)
h = Chem.CanonSmiles('OC(=O)CCC(=O)O')
print('水素付加体（どちらも）', h, Chem.FindMolChiralCenters(Chem.MolFromSmiles(h)))
assert check_equation('C4H2O3 + 2NaOH -> C4H2O4Na2 + H2O') is None
print('ORG-CBX-023 検証OK')
