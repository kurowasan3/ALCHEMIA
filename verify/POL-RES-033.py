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
ITEM = '/home/claude/chemk/outh/items/POL-RES-033/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); check_all_equations(it)
import sympy as sp, random
n, q, rA, rB, x, y = sp.symbols('n q r_A r_B x y', positive=True)
p = rA*x/(rA*x + 1)
mean = sp.simplify((1 - p)*1/(1 - p)**2); assert sp.simplify(mean - (1 + rA*x)) == 0
assert sp.simplify(sp.summation(n*q**(n-1), (n, 1, sp.oo)).args[0][0] - 1/(1-q)**2) == 0 if isinstance(sp.summation(n*q**(n-1), (n, 1, sp.oo)), sp.Piecewise) else True
ML = x*(rA*x + 1)/(rB + x)
assert sp.simplify((1 + rA*x)/(1 + rB/x) - ML) == 0; print('Mayo-Lewis 一致')
v = {rA: 4, rB: sp.Rational(1, 4), x: 1}
pv = p.subs(v); print('p', pv); assert pv == sp.Rational(4, 5)
print('P1', 1 - pv, 'P3', pv**2*(1 - pv)); assert float(pv**2*(1 - pv)) == 0.128
fa = (1 + rA*x)/((1 + rA*x) + (1 + rB/x)); print('A%', fa.subs(v)); assert fa.subs(v) == sp.Rational(4, 5)
xs = sp.solve(sp.Eq(1 + 4*x, 1 + sp.Rational(1, 4)/x), x); print('x', xs); assert xs == [sp.Rational(1, 4)]
print('A in feed', sp.Rational(1, 4)/(1 + sp.Rational(1, 4)))
assert sp.simplify((1 + 2*y)/(1 + sp.Rational(1, 2)/y) - 2*y) == 0
assert sp.simplify(((1 + rA*y)/(1 + rB/y) - rA*y).subs(rB, 1/rA)) == 0
# モンテカルロ（末端モデル）
random.seed(1); pA, pB = 0.8, 0.25/(0.25 + 1)
seq = []; cur = 'A'
for _ in range(2_000_000):
    seq.append(cur)
    cur = ('A' if random.random() < pA else 'B') if cur == 'A' else ('B' if random.random() < pB else 'A')
runs = {'A': [], 'B': []}; L = 1
for i in range(1, len(seq)):
    if seq[i] == seq[i-1]: L += 1
    else: runs[seq[i-1]].append(L); L = 1
ra, rb = runs['A'][1:], runs['B'][1:]
f1 = ra.count(1)/len(ra); f3 = ra.count(3)/len(ra); ma = sum(ra)/len(ra); mb = sum(rb)/len(rb); frac = seq.count('A')/len(seq)
print('MC', round(f1, 3), round(f3, 3), round(ma, 2), round(mb, 3), round(frac, 3))
assert abs(f1 - 0.2) < 0.005 and abs(f3 - 0.128) < 0.005 and abs(ma - 5) < 0.05 and abs(mb - 1.25) < 0.02 and abs(frac - 0.8) < 0.005
print('全検証OK:', ITEM)
