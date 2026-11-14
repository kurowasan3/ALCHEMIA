# 検証スクリプト ORG-STR-037（グループ g01）  python3 -I この.py
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

it = check_item('/home/claude/chem5/out/items/ORG-STR-037/item.json')

import numpy as np, itertools
def rot_group(P):
    G = []
    n = np.cross(P[1] - P[0], P[2] - P[0])
    for perm in itertools.permutations(range(4)):
        Q = P[list(perm)]
        for sgn in (1, -1):
            a = np.cross(P[0], P[1]); planar = np.linalg.norm(a) < 1e-9
            A = np.vstack([P, n if planar else a]); Bm = np.vstack([Q, sgn * (n if planar else np.cross(Q[0], Q[1]))])
            R = np.linalg.lstsq(A, Bm, rcond=None)[0].T
            if np.allclose(R @ A.T, Bm.T, atol=1e-6) and np.allclose(R @ R.T, np.eye(3), atol=1e-6) and np.linalg.det(R) > 0:
                G.append(perm); break
    return G
sq = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float)
th = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
py = np.array([[1, 0, -1], [0, 1, -1], [-1, 0, -1], [0, -1, -1]], float)
def count(P, ms):
    G = rot_group(P); seen = set(); n = 0
    for p in set(itertools.permutations(ms)):
        if p in seen: continue
        n += 1
        for g in G: seen.add(tuple(p[g[i]] for i in range(4)))
    return n
T = {k: [count(P, ms) for P in (sq, th, py)] for k, ms in {'CH2Cl2': 'HHXX', 'CH2ClBr': 'HHXY', 'CHFClBr': 'HXYZ', 'CHCl3': 'HXXX', 'CH3Cl': 'HHHX'}.items()}
print(T)
assert T == {'CH2Cl2': [2, 1, 2], 'CH2ClBr': [2, 1, 3], 'CHFClBr': [3, 2, 6], 'CHCl3': [1, 1, 1], 'CH3Cl': [1, 1, 1]}

print('ALL OK')
