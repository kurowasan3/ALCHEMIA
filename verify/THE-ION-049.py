# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
from scipy.optimize import brentq
ITEM = '/home/claude/chem2/outh/items/THE-ION-049/item.json'
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

Kb, Kw = 2.0e-5, 1.0e-14
def pH(Ka, vb, base='NH3'):
    V = 20.0 + vb; CA = 2.0/V; CB = 0.1*vb/V; KN = Kw/Kb
    if base == 'NH3':
        f = lambda lh: 10**lh + CB*10**lh/(10**lh + KN) - CA*Ka/(Ka + 10**lh) - Kw/10**lh
    else:
        f = lambda lh: 10**lh + CB - CA*Ka/(Ka + 10**lh) - Kw/10**lh
    return -brentq(f, -14, 0)
Ka = 2.0e-4
K = Ka*Kb/Kw; print('(1) K', K); assert abs(K - 4.0e5) < 1e-6
z = 1/(1 + sqrt(K)); print('(2) 未反応', z); assert round(z*100, 2) == 0.16
h = sqrt(Ka*Kw/Kb); print('(2) pH 式', -log10(h), '厳密', pH(Ka, 20.0)); assert round(-log10(h), 2) == 6.50 and abs(pH(Ka, 20.0) - 6.50) < 0.01
# 濃度によらないこと（0.01 mol/L どうしの中和点でも）
print('(3)', pH(Ka, 19.8), pH(Ka, 20.2), 'NaOH', pH(Ka, 20.2, 'NaOH'))
assert round(pH(Ka, 19.8), 1) == 5.7 and round(pH(Ka, 20.2), 1) == 7.3 and round(pH(Ka, 20.2, 'NaOH'), 1) == 10.7
v6 = brentq(lambda v: pH(Ka, v) - 6.0, 19, 21); v7 = brentq(lambda v: pH(Ka, v) - 7.0, 19, 21)
print('(3) 変色範囲', v6, v7); assert round(v6, 1) == 19.9 and round(v7, 1) == 20.1
# (4) 酢酸
Ka2 = 2.0e-5; K2 = Ka2*Kb/Kw; z2 = 1/(1 + sqrt(K2)); print('(4) K', K2, '未反応', z2)
assert abs(K2 - 4.0e4) < 1e-6 and round(z2*100, 2) == 0.50
C = 1.98/39.8; e = 0.020/39.8
y = (-e + sqrt(e*e + 4*C*C/K2))/2; hH = Ka2*(e + y)/C
pe = pH(Ka2, 19.8); print('(4) y', y, '[H+]', hH, 'pH', -log10(hH), '厳密', pe, '近似', -log10(Ka2) + 2)
assert abs(-log10(hH) - pe) < 0.01 and round(hH*1e7, 1) == 2.4 and round(-log10(Ka2) + 2, 2) == 6.70
assert abs(pH(Ka2, 20.0) - 7.00) < 0.01
print('全検証OK:', ITEM)
