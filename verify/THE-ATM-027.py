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
ITEM = '/home/claude/chemk/out/items/THE-ATM-027/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); check_all_equations(it)
import sympy as sp
h, c, eV = sp.Rational(66,10)*sp.Float('1e-34'), sp.Float('3.0e8'), sp.Float('1.6e-19')
# 節の数：2n 個の電子を節 0 から 2 個ずつ詰める
for n in range(1, 11):
    levels = list(range(2*n)); e = 2*n; filled = []
    for k in levels:
        if e <= 0: break
        filled.append(k); e -= 2
    homo, lumo = filled[-1], filled[-1] + 1
    assert (homo, lumo) == (n-1, n)
print('HOMO=n-1, LUMO=n OK; n=3 ->', 2, 3)
Eb = h*c/sp.Float('400e-9')/eV
print('境界 ΔE [eV] =', Eb); assert sig(float(Eb), 3) == 3.09
tab = {4:4.3,5:3.7,6:3.4,7:3.2,8:3.0,9:2.85,10:2.75}
lam = {n: float(h*c/(dE*eV))*1e9 for n, dE in tab.items()}
print({n: round(v,1) for n, v in lam.items()})
nmin = min(n for n, v in lam.items() if v >= 400); assert nmin == 8
assert sig(lam[10], 2) == 450
colors = [(400,435,'黄緑'),(435,480,'黄'),(480,490,'橙'),(490,500,'赤')]
col = [cname for a,b,cname in colors if a <= lam[10] < b]; print('n=10 の色:', col); assert col == ['黄']
assert near(lam[8], 412.5) and near(lam[7], 386.7)
print('全検証OK:', ITEM)
