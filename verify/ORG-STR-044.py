# 検証スクリプト ORG-STR-044（グループ g01）  python3 -I この.py
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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem5/outh/items/ORG-STR-044/item.json')

import numpy as np, itertools
def coords(k, e1, e2):
    th = np.pi / 2 * (k - 1); L = 1.3 * k
    P = [('C%d' % i, np.array([0, 0, 1.3 * i])) for i in range(k + 1)]
    u = np.array([1, 0, 0]); v = np.array([np.cos(th), np.sin(th), 0])
    P += [(e1[0], np.array([0, 0, -0.6]) + u), (e1[1], np.array([0, 0, -0.6]) - u)]
    P += [(e2[0], np.array([0, 0, L + 0.6]) + v), (e2[1], np.array([0, 0, L + 0.6]) - v)]
    return P
def superposable(A, B):
    XA = np.array([x for _, x in A]); order = list(dict.fromkeys(l for l, _ in A))
    idxB = {l: [i for i, (ll, _) in enumerate(B) if ll == l] for l in order}
    idxA = {l: [i for i, (ll, _) in enumerate(A) if ll == l] for l in order}
    for perms in itertools.product(*[list(itertools.permutations(idxB[l])) for l in order]):
        mp = {}
        for l, p in zip(order, perms):
            for ia, ib in zip(idxA[l], p): mp[ia] = ib
        XB = np.array([B[mp[i]][1] for i in range(len(A))]); ca, cb = XA.mean(0), XB.mean(0)
        U, S, Vt = np.linalg.svd((XA - ca).T @ (XB - cb)); d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        if np.allclose((R @ (XA - ca).T).T, XB - cb, atol=1e-6): return True
    return False
mirror = lambda A: [(l, x * np.array([-1, 1, 1])) for l, x in A]
def analyse(k, e1, e2):
    forms = [coords(k, a, b) for a in set(itertools.permutations(e1)) for b in set(itertools.permutations(e2))]
    uniq = []
    for f in forms:
        if not any(superposable(f, u) for u in uniq): uniq.append(f)
    return len(uniq), [not superposable(u, mirror(u)) for u in uniq]
R = {
 'ClHC=C=CHCl': analyse(2, ('Cl', 'H'), ('Cl', 'H')), 'ClHC=C=C=CHCl': analyse(3, ('Cl', 'H'), ('Cl', 'H')),
 'k4': analyse(4, ('Cl', 'H'), ('Cl', 'H')),
 'ア': analyse(2, ('H', 'H'), ('Cl', 'H')), 'イ': analyse(2, ('Me', 'H'), ('Me', 'H')), 'ウ': analyse(2, ('Me', 'Me'), ('Cl', 'H')),
 'エ': analyse(2, ('Me', 'Cl'), ('Me', 'Cl')), 'オ': analyse(3, ('Me', 'H'), ('Me', 'H')),
 'hept-E': analyse(2, ('Me', 'H'), ('E', 'H')), 'hept-Z': analyse(2, ('Me', 'H'), ('Z', 'H')),
 'k1': analyse(1, ('Me', 'H'), ('Et', 'H'))}
for k, v in R.items(): print(k, v)
assert R['ClHC=C=CHCl'] == (2, [True, True]) and R['ClHC=C=C=CHCl'] == (2, [False, False]) and R['k4'] == (2, [True, True])
chiral = [k for k in 'アイウエオ' if any(R[k][1])]
assert chiral == ['イ', 'エ']
assert R['hept-E'] == (2, [True, True]) and R['hept-Z'] == (2, [True, True])   # 2+2=4 種，鏡像の組 2
assert R['k1'] == (2, [False, False])   # 2-ペンテン：シス・トランスとも非キラル

print('ALL OK')
