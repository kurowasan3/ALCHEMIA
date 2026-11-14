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
ITEM = '/home/claude/chemk/out/items/INO-ALK-024/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); assert check_all_equations(it) >= 4
import sympy as sp
H, C, O, Mg, Cl = sp.Rational(10,10), 12, 16, sp.Rational(243,10), sp.Rational(355,10)
MgCO3, MgOH2, W, MgO, MgCl2 = Mg+C+3*O, Mg+2*(O+H), 2*H+O, Mg+O, Mg+2*Cl
M = 3*MgCO3 + MgOH2 + 3*W; print('式量', float(M)); assert near(float(M), 365.2)
assert sig(float(4*MgO/M*100), 3) == 44.1
assert sig(float(sp.Rational(2)/M*3*22400), 3) == 368
# 実験値から逆算
x = sp.Rational(368,22400)/2*100; mg = sp.Rational(441,10)/MgO; y = mg - x
z = (100 - MgCO3*x - MgOH2*y)/W
r = [float(x/y), 1.0, float(z/y)]; print('x:y:z =', r)
assert all(abs(a-b) < 0.05 for a, b in zip(r, [3, 1, 3]))
print('Mg vs CO2 (100g):', float(mg), float(x)); assert mg > x
m = float(sp.Rational(2)/M*4*MgCl2); print('MgCl2 g', m); assert sig(m, 3) == 2.09
print('全検証OK:', ITEM)
