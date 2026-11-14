# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem4/outh/items/THE-ATM-021/item.json'
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
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import numpy as np
from itertools import product
def pat(a, b):
    p = np.poly1d([1])
    for _ in range(a): p = p * np.poly1d([1, 3])
    for _ in range(b): p = p * np.poly1d([1, 1])
    c = p.coeffs[::-1].astype(float)
    return c / c.max() * 100
X = np.array([43, 100, 71, 14]); Zs = np.array([38, 100, 92, 33, 4])
mX = [(a, b) for a, b in product(range(5), range(5)) if a + b == 3 and np.allclose(pat(a, b), X, atol=1.0)]
mZ = [(a, b) for a, b in product(range(5), range(5)) if a + b == 4 and np.allclose(pat(a, b), Zs, atol=1.0)]
print('X', mX, 'Z', mZ); assert mX == [(1, 2)] and mZ == [(2, 2)]
for a, b in product(range(5), range(5)):
    if a + b in (3, 4): print((a, b), np.round(pat(a, b), 1))
def CH(rest, nhal):
    return [(c, rest - 12*c) for c in range(0, rest//12 + 1) if c >= 1 and 0 <= rest - 12*c <= 2*c + 2 - nhal]
assert 206 - 35 - 2*79 == 13 and CH(13, 3) == [(1, 1)]
assert 254 - 2*35 - 2*79 == 26 and CH(26, 4) == [(2, 2)]
r13 = 2*0.011*0.989/0.989**2
print('13C 比', r13); assert abs(r13*100 - 2.2) < 0.05
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from itertools import combinations_with_replacement as cwr
from collections import Counter
allm = Counter({'H': 2, 'Br': 2, 'Cl': 2}); smis = set()
for c1 in cwr(['H', 'Br', 'Cl'], 3):
    c1c = Counter(c1)
    if any(c1c[k] > allm[k] for k in c1c): continue
    c2 = list((allm - c1c).elements())
    f = lambda c: ''.join('(' + t + ')' for t in c if t != 'H')
    smis.add(Chem.MolToSmiles(Chem.MolFromSmiles('C' + f(c1) + 'C' + f(c2))))
stereo = {Chem.MolToSmiles(i) for s in smis for i in EnumerateStereoisomers(Chem.MolFromSmiles(s))}
print(len(smis), len(stereo)); assert len(smis) == 4 and len(stereo) == 6
print('OK')
