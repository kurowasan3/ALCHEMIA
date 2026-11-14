# 検証スクリプト（グループ g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/ORG-HYD-035/item.json'
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
n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b); assert r is None, (b, r); n_eq += 1; print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
from math import log10, floor
def sig(x, n): return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from rdkit import Chem
from sympy import symbols, solve, Rational as R
def ox_sum(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); t = 0
    for a in m.GetAtoms():
        if a.GetSymbol() != 'C': continue
        for b in a.GetBonds():
            o = b.GetOtherAtom(a); bo = int(b.GetBondTypeAsDouble())
            if o.GetSymbol() == 'H': t -= 1
            elif o.GetSymbol() == 'O': t += bo
    return t
def cleave(smi):
    m = Chem.RWMol(Chem.MolFromSmiles(smi))
    b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); m.RemoveBond(i, j)
    prods = []
    for k in (i, j):
        a = m.GetAtomWithIdx(k); h = a.GetTotalNumHs()
        x = m.AddAtom(Chem.Atom(8)); m.AddBond(k, x, Chem.BondType.DOUBLE)
        if h >= 1:
            y = m.AddAtom(Chem.Atom(8)); m.AddBond(k, y, Chem.BondType.SINGLE)
        if h == 2:  # HCOOH → CO2
            pass
    for a in m.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    mm = m.GetMol(); Chem.SanitizeMol(mm)
    frags = [Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mm, asMols=True)]
    return ['O=C=O' if f in ('O=CO',) else f for f in frags]
def nvin(smi):
    m = Chem.MolFromSmiles(smi); b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    return b.GetBeginAtom().GetTotalNumHs() + b.GetEndAtom().GetTotalNumHs()
alk = {'1-ペンテン':'C=CCCC','3-メチル-1-ブテン':'C=CC(C)C','2-ペンテン':'CC=CCC','2-メチル-1-ブテン':'C=C(C)CC','2-メチル-2-ブテン':'CC=C(C)C'}
E = {}
for name, s in alk.items():
    pr = cleave(s); e = sum(ox_sum(p) for p in pr) - ox_sum(s); n = nvin(s)
    assert e == 4 + 2*n, (name, e, n)
    E[name] = (n, e, pr); print(name, 'n =', n, '電子', e, pr)
# (2) 16.0 mL × 0.100 → 1.60 mmol × 5 = 8 e/分子，CO2 なし
cand = [k for k, v in E.items() if v[1] == 8 and 'O=C=O' not in v[2]]
assert cand == ['2-ペンテン'], cand
# (3) 混合物
x = symbols('x'); tot = R(140,1000)/70*R(1,1000)*1000  # mol
tot = R(140,1000)/70
assert tot == R(2,1000)
sol = solve(6*(tot - x) + 8*x - 5*R(1,10)*R(264,10000), x)[0]
assert sol == R(6,10000), sol
print('2-メチル-1-ブテン', float(sol/tot*100), '%  CO2', float(sol*22400), 'mL')
assert sig(float(sol*22400), 3) == 13.4
# (4) 10 e → n=3，水素付加体の一塩素置換体の数
cand4 = [k for k, v in E.items() if v[1] == 10]
assert sorted(cand4) == sorted(['1-ペンテン', '3-メチル-1-ブテン'])
def nclass(smi):
    m = Chem.MolFromSmiles(smi); r = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    return len({r[a.GetIdx()] for a in m.GetAtoms() if a.GetTotalNumHs() > 0})
assert nclass('CCCCC') == 3 and nclass('CCC(C)C') == 4
print('OK')
