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
ITEM = '/home/claude/chemk/outh/items/INO-NIT-030/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); assert check_all_equations(it) >= 8
import sympy as sp
sp_names = ['N2', 'O2', 'NO', 'NO2', 'H2O', 'HNO3']
a = {'N2': -1, 'O2': -1, 'NO': 2}; b = {'NO': -2, 'O2': -1, 'NO2': 2}; c = {'NO2': -3, 'H2O': -1, 'HNO3': 2, 'NO': 1}
net = {k: 2*a.get(k, 0) + 3*b.get(k, 0) + 2*c.get(k, 0) for k in sp_names}; print(net)
assert net == {'N2': -2, 'O2': -5, 'NO': 0, 'NO2': 0, 'H2O': -2, 'HNO3': 4}
N2, O2, Ar = 78.0, 21.0, 1.0
n_used = O2*2/5; rest = N2 - n_used + Ar; print('残り', rest, '減少%', (O2 + n_used)); assert abs(rest - 70.6) < 1e-9
hno3 = O2e = O2*1e-3/22.4*4/5; h = hno3/0.150; import math
pH = -math.log10(h); print('HNO3', hno3, '[H+]', h, 'pH', pH); assert abs(hno3 - 7.5e-4) < 1e-9 and round(pH, 1) == 2.3
need = N2*5/2; print('必要O2', need, '追加', need - O2); assert need - O2 == 174
left = O2 + 200 - need; print('残O2', left); assert left == 26
dec5 = O2 + O2/2; print('(5) 減少', dec5); assert dec5 == 31.5 > O2 + n_used
print('全検証OK:', ITEM)
