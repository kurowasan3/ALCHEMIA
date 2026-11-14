# 検証スクリプト（グループ g03）  python3 -I この.py

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

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_all_ce(it):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    print('検査した反応式:', n)

import sympy as sp
it = json.load(open('/home/claude/chem5/out/items/ORG-CBX-029/item.json', encoding='utf-8'))
check_all_ce(it)
R = sp.Rational
nMn = R(200, 10000) * R(120, 10000)          # 0.0200 mol/L × 12.0 mL
nF = R(5, 2) * nMn
nAcid = R(1, 10) * R(180, 10000)            # 0.100 mol/L × 18.0 mL
nA = nAcid - nF
cF, cA = nF / R(1, 100), nA / R(1, 100)
print('MnO4-', float(nMn), 'HCOOH', float(nF), 'CH3COOH', float(nA), 'c', float(cF), float(cA))
assert nMn == R(24, 100000) and nF == R(6, 10000) and nA == R(12, 10000)
assert float(cF) == 0.06 and float(cA) == 0.12
mAg = 2 * nF * 108
print('Ag', float(mAg)); assert abs(float(mAg) - 0.1296) < 1e-9 and round(float(mAg), 3) == 0.130
M = sp.Matrix([[1, 2], [2, 2]]); print('det', M.det()); assert M.det() != 0
# 電子の数：半反応式
from chemlib import check_equation
for e in ['2MnO4^- + 5HCOOH + 6H+ -> 2Mn^2+ + 5CO2 + 8H2O', 'HCOOH -> CO2 + 2H+ + 2e-', 'MnO4^- + 8H+ + 5e- -> Mn^2+ + 4H2O', '(COOH)2 -> 2CO2 + 2H+ + 2e-']:
    assert check_equation(e) is None, e
print('すべて一致')
