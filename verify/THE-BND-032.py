# 検証スクリプト THE-BND-032  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-BND-032/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
                eqs.add(b)
    return sorted(eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
import itertools
import numpy as np
check_structure()
assert all_equations() == [], all_equations()
s3 = np.sqrt(3)
TBP = {'ax1': (0,0,1), 'ax2': (0,0,-1), 'eq1': (1,0,0), 'eq2': (-.5, s3/2, 0), 'eq3': (-.5, -s3/2, 0)}
OCT = {k: v for k, v in zip('abcdef', [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)])}
def counts(geom, lp):
    ll = lb = 0
    for i, j in itertools.combinations(geom, 2):
        if abs(np.dot(geom[i], geom[j])) < 1e-9:
            a, b = i in lp, j in lp
            if a and b: ll += 1
            elif a or b: lb += 1
    return ll, lb
def best(geom, nlp):
    cands = [(counts(geom, set(c)), c) for c in itertools.combinations(geom, nlp)]
    m = min(c[0] for c in cands)
    return m, [c for k, c in cands if k == m]
def dipole(geom, lp):
    v = sum(np.array(geom[k], float) for k in geom if k not in lp)
    return np.linalg.norm(v)
def bond_angles(geom, lp):
    F = [np.array(geom[k], float) for k in geom if k not in lp]
    return sorted({round(float(np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))))) for a, b in itertools.combinations(F, 2)})
# (1) BrF3
assert counts(TBP, {'ax1','ax2'}) == (0, 6)
assert counts(TBP, {'ax1','eq1'}) == (1, 3)
assert counts(TBP, {'eq1','eq2'}) == (0, 4)
m, sel = best(TBP, 2); print('BrF3', m, sel)
assert m == (0, 4) and all(set(c) <= {'eq1','eq2','eq3'} for c in sel)
assert bond_angles(TBP, set(sel[0])) == [90, 180]       # T字形
# (2) SF4
m, sel = best(TBP, 1); print('SF4', m, sel)
assert m == (0, 2) and all(c[0].startswith('eq') for c in sel)
assert counts(TBP, {'ax1'}) == (0, 3)
assert dipole(TBP, set(sel[0])) > 0.1                    # 極性
# (3) XeF2, XeF4
m, sel = best(TBP, 3); print('XeF2', m, sel)
assert sel == [('eq1','eq2','eq3')] and dipole(TBP, set(sel[0])) < 1e-9
assert bond_angles(TBP, set(sel[0])) == [180]
assert counts(OCT, {'a','c'}) == (1, 6) and counts(OCT, {'a','b'}) == (0, 8)
m, sel = best(OCT, 2); print('XeF4', m, sel)
assert m == (0, 8) and all(np.dot(OCT[c[0]], OCT[c[1]]) < -0.5 for c in sel)
assert dipole(OCT, set(sel[0])) < 1e-9 and bond_angles(OCT, set(sel[0])) == [90, 180]
# (4) IF5
m, sel = best(OCT, 1); assert m == (0, 4) and dipole(OCT, set(sel[0])) > 0.5
# (5) AF3：非共有 k 組で T字形（角 90,180 のみ）かつ v=3+2k<=8
ok = []
for k in range(0, 4):
    n = 3 + k
    if n == 5: geom = TBP
    elif n == 6: geom = OCT
    else: continue   # n=3（正三角形 120°），n=4（三角錐）は 90°/180° を含まない
    m, sel = best(geom, k)
    angs = bond_angles(geom, set(sel[0]))
    print('AF3 k=', k, 'v=', 3+2*k, angs, m)
    if angs == [90, 180] and 3 + 2*k <= 8: ok.append(3 + 2*k)
assert ok == [7], ok
VAL = {'S': 6, 'Br': 7, 'I': 7, 'Xe': 8}
for mol, (A, m_) in {'BrF3': ('Br',3), 'SF4': ('S',4), 'XeF2': ('Xe',2), 'XeF4': ('Xe',4), 'IF5': ('I',5)}.items():
    print(mol, 'n =', m_ + (VAL[A]-m_)//2)
print('全検証OK:', ITEM)
