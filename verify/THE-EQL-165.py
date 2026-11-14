# 検証スクリプト THE-EQL-165   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-EQL-165/item.json'
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
def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def check_structure():
    n_items = len(re.findall(r'\\item\b', re.sub(r'\\begin\{itemize\}.*?\\end\{itemize\}', '', it['problem'], flags=re.S)))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items + 1)), nums
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
check_all_equations()
# ---- 個別の検証（SymPy と数値積分）----
import sympy as sp
k1, km1, k2, M, A, B, X = sp.symbols('k1 km1 k2 M A B X', positive=True)
Xs = sp.solve(sp.Eq(k1*A*B, km1*X + k2*X*M), X)[0]
kp = sp.simplify(k2*Xs*M/(A*B))
assert sp.simplify(kp - k1*k2*M/(km1 + k2*M)) == 0
assert sp.limit(kp/M, M, 0) == k1*k2/km1 and sp.limit(kp, M, sp.oo) == k1
# (3) 2 点からの直線化
d = [(3.0e18, 5.0e-12), (1.5e19, 1.5e-11)]
(x1, y1), (x2, y2) = [(1/m, 1/k) for m, k in d]
slope = (y1 - y2)/(x1 - x2); icpt = y1 - slope*x1
K1 = 1/icpt; Mh = slope*K1
print('k1', K1, 'k-1/k2', Mh); assert sig(K1, 2) == 3.0e-11 and sig(Mh, 2) == 1.5e19
f = lambda m, mh=Mh: K1*m/(m + mh)
g, u = f(2.5e19), f(2.5e18)
print('地表', g, '上空', u, '比', u/g); assert sig(g, 2) == 1.9e-11 and sig(u, 2) == 4.3e-12 and sig(u/g, 2) == 0.23
ar = f(2.5e19, 2*Mh); print('Ar', ar, ar/g); assert sig(ar/g, 2) == 0.73
assert abs(f(1e25, 2*Mh)/f(1e25) - 1) < 1e-5
# 定常状態近似の妥当性：k-1, k2 を具体的に与えて完全な速度式を数値積分し，k' と一致することを確認
import numpy as np
from scipy.integrate import solve_ivp
km1v = 1.0e9; k2v = km1v/Mh; Mv = 2.5e19
def rhs(t, y):
    oh, no2, x, h = y
    r1 = K1*oh*no2; rm = km1v*x; r2 = k2v*x*Mv
    return [-r1 + rm, -r1 + rm, r1 - rm - r2, r2]
y0 = [1.0e7, 1.0e11, 0, 0]
sol = solve_ivp(rhs, [0, 2.0], y0, method='LSODA', rtol=1e-8, atol=1e-6, dense_output=True)
t = 1.0; oh = sol.sol(t)[0]
kobs = -np.log(oh/1.0e7)/(1.0e11*t)
print('数値積分の見かけの k', kobs, '式', f(Mv)); assert abs(kobs/f(Mv) - 1) < 0.01
print('すべて OK')
