# 検証スクリプト THE-BND-035  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-BND-035/item.json'
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
from sympy import sqrt, Rational as R, nsimplify
check_structure()
assert all_equations() == []
ideal = 2*sqrt(R(2,3)); assert ideal == sqrt(R(8,3)); assert sig(float(ideal), 3) == 1.63
assert sig(2*1.414/1.73, 3) == 1.63
# 実際に hcp 格子を作って最近接距離を数える
def neigh(a, c):
    A = np.array([a, 0, 0]); B = np.array([a/2, a*np.sqrt(3)/2, 0]); C = np.array([0, 0, c])
    basis = [np.zeros(3), (A + B)/3 + C/2]   # 60°基底で (1/3,1/3,1/2)
    pts = [i*A + j*B + k*C + b for i, j, k in itertools.product(range(-2, 3), repeat=3) for b in basis]
    d = sorted(round(float(np.linalg.norm(p)), 4) for p in pts if np.linalg.norm(p) > 1e-9)
    out = {}
    for x in d[:12]: out[x] = out.get(x, 0) + 1
    return out
data = {'Mg': (0.321, 0.521), 'Zn': (0.266, 0.495), 'Be': (0.229, 0.358)}
res = {}
for m, (a, c) in data.items():
    n = neigh(a, c); dlay = np.sqrt(a*a/3 + c*c/4); print(m, n, round(c/a, 4), round(dlay, 4))
    assert n.get(round(a, 4)) == 6 and n.get(round(float(dlay), 4)) == 6
    r = min(a, dlay)/2
    phi_exact = 6*4/3*np.pi*r**3 / (3*np.sqrt(3)/2*a*a*c)
    phi_given = 6*4/3*3.14*r**3 / (3*1.73/2*a*a*c)
    res[m] = (r, phi_exact, phi_given, dlay); print('  r', r, 'phi', phi_exact, phi_given)
assert sig(res['Zn'][0], 3) == 0.133 and sig(res['Be'][0], 3) == 0.111
assert data['Zn'][1]/data['Zn'][0] > float(ideal) > data['Be'][1]/data['Be'][0]
for m, p in (('Zn', 65), ('Be', 71), ('Mg', 74)):
    assert round(res[m][1]*100) == p and round(res[m][2]*100) == p, m
assert sig(res['Mg'][3], 3) == 0.320
a, c = data['Zn']
V = 3*1.73/2*(a*1e-7)**2*(c*1e-7); rho = 6*65.4/(6.0e23*V); print('V', V, 'rho', rho)
assert sig(V, 3) == 9.09e-23 and sig(rho, 2) == 7.2
rho_exact = 6*65.4/(6.02e23*3*np.sqrt(3)/2*(0.2665e-7)**2*0.4947e-7); print('rho(文献格子定数)', rho_exact)
print('全検証OK:', ITEM)
