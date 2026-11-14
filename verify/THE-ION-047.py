# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
from scipy.optimize import brentq
from sympy import symbols, solve, Rational
ITEM = '/home/claude/chem2/outh/items/THE-ION-047/item.json'
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
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and ' + ' in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

KB, Kb, Kw = 6.0e-10, 2.0e-5, 1.0e-14
K = KB*Kb/Kw; print('(1) K =', K); assert abs(K - 1.2) < 1e-12
x = symbols('x', positive=True)
n = 0.0400
res = {}
for C in (0.400, 0.0400):
    r = [s for s in solve(x**2 - 1.2*(n - x)*(C - x), x) if s < n][0]
    res[C] = (float(r), float((n - r)/n)); print('(2) C=', C, 'x=', r, '分子の割合', (n - r)/n)
assert round(res[0.400][1]*100, 1) == 7.8 and round(res[0.0400][1]*100) == 48
assert round(0.0031/0.0400*100, 1) == 7.8 and round(0.0191/0.0400*100) == 48

# 厳密（H+・OH- を含む）: 吸収液
def absorb(C, n, V_add=0.0, cHCl=0.100):
    V = 50.0 + V_add
    N = n*50.0/V; B = C*50.0/V; Cl = cHCl*V_add/V
    KN = Kw/Kb
    def f(lh):
        h = 10**lh
        NH4 = N*h/(h + KN); BOH4 = B*KB/(KB + h)
        return h + NH4 - BOH4 - Cl - Kw/h
    lh = brentq(f, -14, 0); h = 10**lh
    return -lh, N*KN/(h + KN)/N
p0, free = absorb(0.400, 0.0400); print('厳密 pH', p0, '分子の割合', free)
assert abs(free - res[0.400][1]) < 0.002 and round(p0, 1) == 8.2

# (3)(4) 滴定曲線
V_eq = 2.00/0.100; print('(3)', V_eq)
p_eq, _ = absorb(0.400, 0.0400, 20.0)
h4 = sqrt(KB*20/70 + (Kw/Kb)*2/70); print('(4) 式', -log10(h4), '厳密', p_eq)
assert round(-log10(h4), 2) == 4.87 and abs(p_eq - 4.87) < 0.01
assert abs(sqrt(1.86) - 1.36) < 0.005
# 当量点の前後で pH が大きく変わること
print('19.9 mL', absorb(0.4, 0.04, 19.9)[0], '20.1 mL', absorb(0.4, 0.04, 20.1)[0])
# (5) pH 4.0 になる滴下量
v4 = brentq(lambda v: absorb(0.4, 0.04, v)[0] - 4.0, 20.0, 21.0)
print('(5) pH4.0 の滴下量', v4, '過剰', v4 - 20.0, '近似', 1.0e-4*0.0700/0.100*1000)
assert abs((v4 - 20.0) - 0.070) < 0.01 and round(1.0e-4*0.0700/0.100*1000, 3) == 0.070
print('全検証OK:', ITEM)
