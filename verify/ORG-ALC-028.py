# 検証スクリプト（グループ g03）  python3 -I この.py

import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_all_ce(it):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    print('検査した反応式:', n)

from rdkit import Chem
from itertools import combinations
import sympy as sp
it = json.load(open('/home/claude/chem5/outh/items/ORG-ALC-028/item.json', encoding='utf-8'))
check_all_ce(it)

def can(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s))
names = {'OCCCCC':'1-ペンタノール','CC(O)CCC':'2-ペンタノール','CCC(O)CC':'3-ペンタノール','OCC(C)CC':'2-メチル-1-ブタノール',
 'CC(C)(O)CC':'2-メチル-2-ブタノール','CC(O)C(C)C':'3-メチル-2-ブタノール','OCCC(C)C':'3-メチル-1-ブタノール','OCC(C)(C)C':'2,2-ジメチル-1-プロパノール',
 'CC(O)CC':'2-ブタノール'}
names = {can(k): v for k, v in names.items()}

def betaH(s):
    """[(アルケンSMILES, 置換度k, シストランスの有無, Hの数)]"""
    m = Chem.MolFromSmiles(s)
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]; c = o.GetNeighbors()[0]
    out = []
    for b in c.GetNeighbors():
        if b.GetSymbol() != 'C' or b.GetTotalNumHs() == 0: continue
        rw = Chem.RWMol(m); ci, bi, oi = c.GetIdx(), b.GetIdx(), o.GetIdx()
        rw.RemoveAtom(oi); ci -= (oi < ci); bi -= (oi < bi)
        rw.GetBondBetweenAtoms(ci, bi).SetBondType(Chem.BondType.DOUBLE)
        p = rw.GetMol(); Chem.SanitizeMol(p)
        k = sum(1 for a in (ci, bi) for n in p.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() not in (ci, bi))
        ct = all(len({Chem.CanonicalRankAtoms(p, breakTies=False)[n.GetIdx()] for n in p.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() not in (ci, bi)} | ({'H'} if p.GetAtomWithIdx(a).GetTotalNumHs() else set())) == 2 for a in (ci, bi))
        out.append((Chem.MolToSmiles(p), k, ct, b.GetTotalNumHs()))
    return out

f1, f2, f3, r = sp.symbols('f1 f2 f3 r', positive=True)
F = {1: f1, 2: f2, 3: f3}
def comp(s, vals):
    w = {}
    for alk, k, ct, nh in betaH(s):
        wt = nh * F[k]
        if ct:
            w[('trans', alk)] = w.get(('trans', alk), 0) + wt * r / (r + 1)
            w[('cis', alk)] = w.get(('cis', alk), 0) + wt / (r + 1)
        else:
            w[alk] = w.get(alk, 0) + wt
    tot = sum(w.values())
    return {kk: sp.nsimplify(sp.simplify((v / tot).subs(vals))) for kk, v in w.items()} if w else {}

# 実験 1・2 から解く
b2 = betaH(can('CC(O)CC')); print('2-ブタノール', b2)
sol = sp.solve([sp.Eq(2*f2/(3*f1), 4), sp.Eq(r, 3), sp.Eq(2*f3/(6*f2), sp.Rational(7, 3))], [f2, f3, r], dict=True)[0]
vals = {f1: 1, **{k: v.subs(f1, 1) for k, v in sol.items()}}
print('f2/f1 =', vals[f2], 'f3/f1 =', vals[f3], 'r =', vals[r])
assert vals[f2] == 6 and vals[f3] == 42 and vals[r] == 3
c2b = comp(can('CC(O)CC'), vals); print(c2b)
assert sorted(c2b.values()) == [sp.Rational(1,5), sp.Rational(1,5), sp.Rational(3,5)]
c22 = comp(can('CC(C)(O)CC'), vals); print(c22)
# (2)
c = comp(can('CC(O)C(C)C'), vals); print('3-メチル-2-ブタノール', c)
assert c[can('C=CC(C)C')] == sp.Rational(1, 15) and c[can('CC=C(C)C')] == sp.Rational(14, 15)
print('6.7% / 93%:', float(c[can('C=CC(C)C')])*100, float(c[can('CC=C(C)C')])*100)
# (3) 8 種の 2 種組み合わせを総当たり
C5 = [s for s in names if names[s] != '2-ブタノール']
comps = {s: comp(s, vals) for s in C5}
for s in C5: print(names[s], {k: float(v) for k, v in comps[s].items()})
obs = {can('C=C(C)CC'): 0.100, can('CC=C(C)C'): 0.856, can('C=CC(C)C'): 0.044}
hits = []
import numpy as np
for a, b in combinations(C5, 2):
    keys = set(comps[a]) | set(comps[b]) | set(obs)
    best = None
    for x in np.linspace(0, 1, 100001):
        na, nb = x, 1 - x
        # 脱水しないアルコールはアルケンを生じない
        tot = na * (1 if comps[a] else 0) + nb * (1 if comps[b] else 0)
        if tot == 0: continue
        pred = {k: (na * float(comps[a].get(k, 0)) + nb * float(comps[b].get(k, 0))) / tot for k in keys}
        err = max(abs(pred[k] - obs.get(k, 0)) for k in keys)
        if best is None or err < best[0]: best = (err, x)
    if best[0] < 0.0006:
        hits.append((names[a], names[b], best))
print(hits)
assert len(hits) == 1 and {hits[0][0], hits[0][1]} == {'2-メチル-2-ブタノール', '3-メチル-2-ブタノール'}
x = hits[0][2][1] if hits[0][0] == '2-メチル-2-ブタノール' else 1 - hits[0][2][1]
assert abs(x - 1/3) < 1e-3
print('物質量比 1:2 を確認．すべて一致')
