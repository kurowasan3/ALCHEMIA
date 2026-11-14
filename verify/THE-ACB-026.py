# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
from scipy.optimize import brentq
ITEM = '/home/claude/chem2/outh/items/THE-ACB-026/item.json'
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

K1, K2, Kw = 4.0e-7, 5.0e-11, 1.0e-14
def pH(nC, nNa, nCl, V):
    C, Na, Cl = nC / V, nNa / V, nCl / V
    def f(lh):
        h = 10**lh; D = h*h + K1*h + K1*K2
        return h + Na - Cl - C*K1*h/D - 2*C*K1*K2/D - Kw/h
    return -brentq(f, -14, 0)
# (1)(2) CO2 がとどまる
p = {v: pH(1.0, 2.0, 0.1*v, 20 + v) for v in (19.8, 20.0, 20.2)}
print('CO2残', p)
assert round(-log10(sqrt(K1*0.025)), 1) == 4.0 and abs(p[20.0] - 4.0) < 0.01
assert round(6.40 - 1.70, 1) == 4.7 and round(p[19.8], 1) == 4.7
assert round(p[20.2], 1) == 3.3
# 近似の確認
h198 = K1*49; print('H+/HCO3-', h198 / (0.02/39.8)); assert h198 / (0.02/39.8) < 0.05
print('H2CO3 由来/過剰H+', K1*0.0249/5e-4 / 5e-4)
# (3) 煮沸後
C = 0.020/40.0
hf = sqrt(K1*(K2*C + Kw)/(K1 + C)); pf = -log10(hf)
pe = pH(0.020, 2.0, 1.98, 40.0)
print('式', pf, '厳密', pe, '近似', (-log10(K1) - log10(K2))/2)
assert round(pf, 2) == 8.28 and abs(pe - pf) < 0.01 and round((-log10(K1) - log10(K2))/2, 2) == 8.35
assert abs(0.5*(17 - 0.45) - pf) < 0.01
# (4)
p200 = pH(0.0, 2.0, 2.00, 40.0); p202 = pH(0.0, 2.0, 2.02, 40.0)
print('煮沸後 20.0', p200, '20.2', p202); assert round(p202, 1) == 3.3 and abs(p200 - 7) < 0.01
# (5) 煮沸で HCO3- の一部が CO3^2- + CO2 となっても，終点（H+ 受け取り総量）は Na+ で決まる
#   例：残り 0.020 mmol の半分が反応したとき，0.005 mmol CO3^2- + 0.010 mmol HCO3- → 必要な H+ = 2*0.005+0.010 = 0.020 mmol
assert abs(2*0.005 + 0.010 - 0.020) < 1e-12
print('全検証OK:', ITEM)
