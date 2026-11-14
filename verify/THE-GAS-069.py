# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-GAS-069/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
NA = 6.02e23; Vm = 22.4e3  # cm3/mol
n_d = 4.0e12; H = 25e5     # cm
N = n_d*H; n = N/NA
assert sig(N, 2) == 1.0e19 and sig(n, 2) == 1.7e-5
th = n*Vm*10   # mm
print('厚さ mm', th); assert sig(th, 2) == 3.7
c = n_d/NA*1000; bc = H*c
print('c', c, 'bc', bc); assert sig(c, 2) == 6.6e-9 and sig(bc, 2) == 1.7e-2
bc_stp = (1/22.4)*(n*Vm)
assert abs(bc_stp - bc)/bc < 1e-9
A255 = 2950*bc; A310 = 24*bc
print(A255, A310); assert round(A255) == 49 and round(A310, 2) == 0.40
assert round(1 - 10**-A310, 2) == 0.60 and abs(10**-0.40 - 0.398) < 5e-4 and abs(10**-0.36 - 0.437) < 5e-4
A310b = 0.9*A310; A255b = 0.9*A255
assert round(A255b) == 44 and round(A310b, 2) == 0.36
ratio = 10**-A310b/10**-A310; print('増加', ratio); assert round(ratio, 2) == 1.10
# ε と吸収断面積の対応（実測の目安：255 nm で σ≈1.1e-17 cm2 → ε=σNA/2303≈2.9e3）
print('ε(255) from σ', 1.1e-17*NA/2303)
print('OK')
