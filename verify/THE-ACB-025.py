# 検証スクリプト THE-ACB-025  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-ACB-025/item.json'
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
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    return re.sub(r'\{(\d*[+-])\}', r'\1', eq)
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'<=>|->|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                e = check_equation(clean_ce(b)); n += 1
                assert e is None, (b, e)
    print('反応式', n, '本 つり合いOK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

from sympy import symbols, nsolve, log as slog
from math import log10, sqrt
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
rule = lambda m: 8 - 5*m
assert (rule(0), rule(2), rule(3)) == (8, -2, -7)
pred = [rule(1) + 5*i for i in range(3)]; assert pred == [3, 8, 13]
assert all(abs(p - e) <= 1.0 for p, e in zip(pred, [2.1, 7.2, 12.3]))
# 構造の列挙: H3POx で O を「OH」「=O」に分け，残りの H を P に結合させる
def cands(nO, nH, n_eq):
    out = []
    for nOH in range(0, nO+1):
        m = nO - nOH; nPH = nH - nOH
        if nPH < 0: continue
        # P の結合数: 単結合 OH, P-H, 二重結合 =O で原子価 5（4配位）
        if nOH + nPH + m == 4 and nOH + nPH + 2*m == 5 and nOH == n_eq: out.append((m, nOH, nPH))
    return out
c3 = cands(3, 3, 2); c2 = cands(2, 3, 1); print(c3, c2)
assert c3 == [(1, 2, 1)] and c2 == [(1, 1, 2)]
assert abs(rule(1) - 1.3) < 2 and abs(rule(0) - 1.3) > 6
ps = Chem.SmilesParserParams(); ps.removeHs = False
assert CalcMolFormula(Chem.MolFromSmiles('P(=O)(O)(O)[H]', ps)) == 'H3O3P'
assert CalcMolFormula(Chem.MolFromSmiles('P(=O)(O)([H])[H]', ps)) == 'H3O2P'
ox = lambda nH, nO: 2*nO - nH
assert (ox(3,4), ox(3,3), ox(3,2)) == (5, 3, 1)
K1, K2, c, Kw = 5.0e-2, 2.0e-7, 0.0500, 1e-14
h = symbols('h', positive=True); D = h**2 + K1*h + K1*K2
ex = nsolve(h + c - c*(K1*h + 2*K1*K2)/D - Kw/h, h, 7e-5); pH_ex = float(-slog(ex, 10)); print('厳密', pH_ex)
pH_ap = -log10(sqrt(K1*K2*c/(K1+c))); print('近似', pH_ap)
assert round(pH_ex, 2) == 4.15 and round(pH_ap, 2) == 4.15
assert abs((1.30+6.70)/2 - 4.00) < 1e-9
assert abs(-0.5*(0.70-9) - 4.15) < 1e-9
r = float(ex)/K1; print('H2A/HA-', r); assert r < 2e-3

print('全検証OK:', ITEM)
