# 検証スクリプト THE-SOL-111   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-SOL-111/item.json'
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
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import sympy as sp
p0, V0 = sp.Rational(28,10), 1
P1 = p0*V0/(V0-sp.Rational(1,2)); W1 = P1*sp.Rational(1,2)
print('(1)', P1, W1); assert float(P1) == 5.6 and float(W1) == 2.8
v = sp.symbols('v')
W2 = p0*v/(1-v) + P1*(sp.Rational(1,2)-v)
print('(2) 等分', float(W2.subs(v, sp.Rational(1,4)))); assert sig(float(W2.subs(v, sp.Rational(1,4))),3) == 2.33
vs = [s for s in sp.solve(sp.diff(W2, v), v) if 0 < s < sp.Rational(1,2)][0]
print('最適 v1', float(vs), float(W2.subs(v, vs))); assert sig(float(vs),3) == 0.293 and sig(float(W2.subs(v, vs)),3) == 2.32
# 数値で最小であることの確認
import numpy as np
g = sp.lambdify(v, W2); xs = np.linspace(0.01, 0.49, 4801); assert abs(xs[np.argmin(g(xs))] - float(vs)) < 1e-3
Winf = p0*sp.log(2); print('(3)', float(Winf), float(W1/Winf)); assert sig(float(Winf),3) == 1.94 and sig(float(W1/Winf),3) == 1.44
# 無限段の和の極限（段数 N の階段近似が W∞ に収束）
N = 2000; edges = np.linspace(0, 0.5, N+1); stair = sum(2.8/(1-edges[i+1])*(edges[i+1]-edges[i]) for i in range(N))
print('階段近似', stair); assert abs(stair - float(Winf)) < 1e-2
r = sp.symbols('r'); w = p0/(r*(1-r)); rs = sp.solve(sp.diff(w, r), r)
print('(4)', rs, float(w.subs(r, rs[0]))); assert rs == [sp.Rational(1,2)] and float(w.subs(r, rs[0])) == 11.2
assert sp.simplify(sp.diff(p0/(1-r), r) - p0/(1-r)**2) == 0  # 単調増加

print('検証OK')
