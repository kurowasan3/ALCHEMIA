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

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import re as _re
it = json.load(open('/home/claude/chem5/out/items/THE-RDX-031/item.json', encoding='utf-8'))
check_all_ce(it)
def oxsum(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); tot = 0
    for a in m.GetAtoms():
        if a.GetSymbol() != 'C': continue
        for b in a.GetBonds():
            o = b.GetOtherAtom(a); order = int(b.GetBondTypeAsDouble())
            if o.GetSymbol() == 'H': tot -= 1
            elif o.GetSymbol() in ('O', 'Br', 'Cl'): tot += order
    return tot
def byformula(smi):
    f = CalcMolFormula(Chem.MolFromSmiles(smi)); cnt = dict((k, int(v) if v else 1) for k, v in _re.findall(r'([A-Z][a-z]?)(\d*)', f))
    return 2*cnt.get('O', 0) + cnt.get('Br', 0) + cnt.get('Cl', 0) - cnt.get('H', 0)
mols = {'A': 'CC(O)C', 'B': 'C=CC', 'C': 'CC(=O)C', 'D': 'BrCC(Br)C', 'E': 'CCC', 'F': 'CC(Cl)C', '1-プロパノール': 'CCCO', 'エチルメチルエーテル': 'CCOC'}
res = {k: oxsum(v) for k, v in mols.items()}; print(res)
assert res == {'A': -6, 'B': -6, 'C': -4, 'D': -4, 'E': -8, 'F': -6, '1-プロパノール': -6, 'エチルメチルエーテル': -6}
assert all(oxsum(v) == byformula(v) for v in mols.values())
judge = {'ア': ('B', 'A'), 'イ': ('A', 'C'), 'ウ': ('C', 'A'), 'エ': ('B', 'D'), 'オ': ('B', 'E'), 'カ': ('E', 'F'), 'キ': ('A', 'B')}
out = {k: ('酸化' if res[b] > res[a] else '還元' if res[b] < res[a] else 'なし') for k, (a, b) in judge.items()}
print(out)
assert out == {'ア': 'なし', 'イ': '酸化', 'ウ': '還元', 'エ': '酸化', 'オ': '還元', 'カ': '酸化', 'キ': 'なし'}
e = 3*4 - res['A']; print('電子', e, 'O2', e/4); assert e == 18 and e/4 == 4.5
assert check_equation('2C3H8O + 9O2 -> 6CO2 + 8H2O') is None
print('すべて一致')
