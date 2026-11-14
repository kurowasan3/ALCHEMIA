# 検証スクリプト THE-EQL-164   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-EQL-164/item.json'
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
# ---- 個別の検証（SymPy）----
import sympy as sp
x, r, P = sp.symbols('x r P', positive=True)
n = {'SiHCl3': 1 - x, 'H2': r - x, 'HCl': 3 * x}; nt = sum(n.values())
p = {k: v / nt * P for k, v in n.items()}
Kexpr = sp.simplify(p['HCl']**3 / (p['SiHCl3'] * p['H2']))
assert sp.simplify(Kexpr - 27 * x**3 * P / ((1 - x) * (r - x) * (1 + r + x))) == 0
K = float(Kexpr.subs({x: 0.30, r: 10, P: 1.0e5})); print('Kp', K); assert sig(K, 2) == 950
f = lambda xx, rr, PP: float(Kexpr.subs({x: xx, r: rr, P: PP}))
print('f(0.36)', f(0.36, 10, 5e4), 'f(0.37)', f(0.37, 10, 5e4))
assert sig(f(0.36, 10, 5e4), 2) == 900 and sig(f(0.37, 10, 5e4), 2) == 990
for Kv in (K, 950.0):
    x3 = float(sp.nsolve(Kexpr.subs({r: 10, P: 5e4}) - Kv, x, 0.35)); print('x(P=5e4)', x3); assert sig(x3, 2) == 0.37
    xa = float(sp.nsolve(27 * x**3 * 5e4 / ((1 - x) * 10 * 11) - Kv, x, 0.35)); print('近似', xa); assert sig(xa, 2) == 0.37
    x4 = float(sp.nsolve(Kexpr.subs({r: 20, P: 1e5}) - Kv, x, 0.4)); print('x(r=20)', x4); assert sig(x4, 2) == 0.44
assert sig(0.30 / 10, 2) == 0.030 and sig(0.44 / 20, 2) == 0.022
Q1 = (5.0e3)**3 / (1.0e4 * 8.5e4); Q2 = (1.8e4)**3 / (2.0e3 * 8.0e4)
print('Q1', Q1, 'Q2', Q2); assert Q1 < 950 < Q2 and sig(Q1, 2) == 150 and sig(Q2, 2) == 36000
# 分圧の和が全圧に近い（1 気圧程度の反応器）ことの確認
print('(I) の全圧', 1.0e4 + 8.5e4 + 5.0e3, '(II) の全圧', 2.0e3 + 8.0e4 + 1.8e4)
print('すべて OK')
