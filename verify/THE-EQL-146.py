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
ITEM = '/home/claude/chemk/out/items/THE-EQL-146/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); assert check_all_equations(it) >= 1
data = [(2.0e-3, 1.0e-3, 1.2e-6), (4.0e-3, 3.0e-3, 1.44e-5), (6.0e-3, 2.0e-3, 2.16e-5)]
sols = [(x, y) for x in range(5) for y in range(5) if abs(2**x*3**y - 12) < 1e-9 and abs(3**x*2**y - 18) < 1e-9]
print('解', sols); assert sols == [(2, 1)]
ks = [v/(a**2*b) for a, b, v in data]; print('k', ks); assert all(near(k, 300) for k in ks)
k = 300.0; no0, h0 = 1.0e-3, 5.0e-2
est = 2*k*no0**2*h0*2.0/no0*100; print('見積もり %', est); assert sig(est, 2) == 6.0
# 数値積分（RK4）
def f(s):
    no, h2 = s; r = k*no*no*h2; return (-2*r, -2*r)
s = (no0, h0); dt = 1e-4
for _ in range(int(2.0/dt)):
    k1 = f(s); k2 = f(tuple(a + dt/2*b for a, b in zip(s, k1))); k3 = f(tuple(a + dt/2*b for a, b in zip(s, k2))); k4 = f(tuple(a + dt*b for a, b in zip(s, k3)))
    s = tuple(a + dt/6*(b + 2*c + 2*d + e) for a, b, c, d, e in zip(s, k1, k2, k3, k4))
real = (no0 - s[0])/no0*100; print('実際 %', real); assert real < est and abs(real - 5.66) < 0.02
print('全検証OK:', ITEM)
