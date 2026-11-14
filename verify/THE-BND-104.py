# 検証スクリプト（グループ g02）  python3 -I この.py
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

def check_all_equations(it, extra=()):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|…', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]', b): continue
            if not (' + ' in b or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', b)): continue
            e = clean_ce(b)
            r = check_equation(e)
            assert r is None, (b, r)
            print('式OK:', b); n += 1
    for e in extra:
        r = check_equation(e); assert r is None, (e, r); print('式OK(追加):', e); n += 1
    return n

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def check_parts(it):
    k = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, k + 1)], (k, nums)
    print('小問数OK:', k)
ITEM = '/home/claude/chemk/outh/items/THE-BND-104/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); check_all_equations(it)
import numpy as np
from fractions import Fraction as Fr
a, c = 0.3029, 0.3228
A = np.array([[a, 0, 0], [-a/2, a*np.sqrt(3)/2, 0], [0, 0, c]])   # 格子ベクトル（行）
atoms = [('Ti', (0, 0, 0)), ('B', (1/3, 2/3, 1/2)), ('B', (2/3, 1/3, 1/2))]
def neighbors(i):
    _, fi = atoms[i]; ri = np.array(fi) @ A; out = []
    for j, (sj, fj) in enumerate(atoms):
        for n1 in range(-2, 3):
            for n2 in range(-2, 3):
                for n3 in range(-2, 3):
                    r = (np.array(fj) + [n1, n2, n3]) @ A
                    dd = np.linalg.norm(r - ri)
                    if dd > 1e-9: out.append((round(dd, 4), sj))
    return sorted(out)
nb = neighbors(1)
bb = sorted(d for d, s in nb if s == 'B')[:3]; bt = [d for d, s in nb if s == 'Ti']; bt = [d for d in bt if abs(d - min(bt)) < 1e-3]
print('B-B', bb, 'B-Ti', len(bt), bt[0]); assert len(set(bb)) == 1 and sig(bb[0], 3) == 0.175 and len(bt) == 6 and sig(bt[0], 3) == 0.238
nt = neighbors(0); tb = [d for d, s in nt if s == 'B']; tb = [d for d in tb if abs(d - min(tb)) < 1e-3]
print('Ti-B', len(tb)); assert len(tb) == 12
assert 2*Fr(1,12) + 2*Fr(1,6) == Fr(1,2)
V = np.sqrt(3)/2*a*a*c*1e-21; rho = (47.9 + 2*10.8)/6.02e23/V
V2 = 1.732/2*a*a*c*1e-21; rho2 = (47.9 + 2*10.8)/6.02e23/V2
print('V', V, 'rho', rho, rho2); assert sig(rho, 3) == 4.50 and sig(rho2, 3) == 4.50
print('全検証OK:', ITEM)
