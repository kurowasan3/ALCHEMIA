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
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it = json.load(open('/home/claude/chem5/outh/items/ORG-CBX-036/item.json', encoding='utf-8'))
check_all_ce(it)
# 一般式 R はエタノールに置きかえて確認
assert check_equation('C2H5OH + (CH3CO)2O -> CH3COOC2H5 + CH3COOH') is None
assert check_equation('C3H5(OH)3 + 3(CH3CO)2O -> C3H5(OCOCH3)3 + 3CH3COOH') is None
R = sp.Rational
c = 1; V0 = R(5000, 100)  # mL
nAc2O = c * V0 / 1000 / 2; print('無水酢酸', nAc2O); assert nAc2O == R(25, 1000)
nOH_A = c * (V0 - R(3500, 100)) / 1000; print('A の OH', nOH_A); assert nOH_A == R(15, 1000) and nOH_A < nAc2O
sols = []
for m in range(1, 13):
    M = R(460, 1000) * m / nOH_A
    if M > 200: break
    n = (M - 2 - 16*m) / 14
    print('m', m, 'M', M, 'n', n)
    if n.q == 1 and n >= m: sols.append((int(n), m))
assert sols == [(3, 3)], sols
assert CalcMolFormula(Chem.MolFromSmiles('OCC(O)CO')) == 'C3H8O3'
# B
nCOOH = c * R(500, 100) / 1000
nOH_B = c * (V0 - V0) / 1000 + nCOOH
M_B = R(450, 1000) / nCOOH
print('B: COOH', nCOOH, 'OH', nOH_B, 'M', M_B); assert nOH_B == nCOOH and M_B == 90
# 分子量 90 で O>=3 の C,H,O の式（H は 2C+2 以下，偶数）
cands = [(x, y, z) for z in range(3, 6) for x in range(1, 8) for y in range(0, 2*x+3) if 12*x + y + 16*z == 90 and y % 2 == 0]
print('C,H,O の候補', cands)
# COOH と OH が 1 個ずつ（O が 3 個必要）→ O4 の C2H2O4 は H が COOH 2 個分しかない（シュウ酸）
assert (3, 6, 3) in cands
# C3H6O3 の鎖式で COOH 1 個・OH 1 個をもつもの
c3 = ['CC(O)C(=O)O', 'OCCC(=O)O']
for s in c3:
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C3H6O3'
chir = [s for s in c3 if Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False)]
assert chir == ['CC(O)C(=O)O']
print('すべて一致')
