# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/out/items/ORG-STR-019/item.json'
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
R = sp.Rational
mC = R(968, 10) * 12 / 44; mH = R(342, 10) * 2 / 18; mO = R(510, 10) - mC - mH
r = [mC / 12, mH / 1, mO / 16]
base = r[2] / 13
ratio = [x / base for x in r]
print('C:H:O =', ratio)
assert ratio == [22, 38, 13]
M = 22 * 12 + 38 + 13 * 16
assert M == 510 and 2 * M > 600
k, m = sp.symbols('k m', positive=True, integer=True)
sol = sp.solve([12 + k * m - 22, 22 + 2 * k * m - 2 * k - 38, 11 + 2 * k - k - 13], [k, m], dict=True)
print(sol); assert sol == [{k: 2, m: 5}]
acids = ['CCCCC(=O)O', 'CCC(C)C(=O)O', 'CC(C)CC(=O)O', 'CC(C)(C)C(=O)O']
assert len({Chem.CanonSmiles(a) for a in acids}) == 4
for a in acids:
    mm = Chem.MolFromSmiles(a)
    assert CalcMolFormula(mm) == 'C5H10O2'
    print(a, Chem.FindMolChiralCenters(mm, includeUnassigned=True))
chir = [a for a in acids if Chem.FindMolChiralCenters(Chem.MolFromSmiles(a), includeUnassigned=True)]
assert chir == ['CCC(C)C(=O)O']
# 例：スクロースの 6 位と 6' 位を 2-メチルブタン酸でエステル化
suc = 'OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O'
assert CalcMolFormula(Chem.MolFromSmiles(suc)) == 'C12H22O11'
est = 'CCC(C)C(=O)OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](COC(=O)C(C)CC)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O'
assert CalcMolFormula(Chem.MolFromSmiles(est)) == 'C22H38O13'
print('ORG-STR-019 検証OK')
