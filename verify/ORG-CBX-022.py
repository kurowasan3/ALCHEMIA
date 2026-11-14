# 検証スクリプト ORG-CBX-022  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/ORG-CBX-022/item.json'
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

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
pKa = {'FCH2COOH': (2.59, 'FCC(=O)O', 'C2H3FO2'), 'ClCH2COOH': (2.87, 'ClCC(=O)O', 'C2H3ClO2'),
       'BrCH2COOH': (2.90, 'BrCC(=O)O', 'C2H3BrO2'), 'ICH2COOH': (3.18, 'ICC(=O)O', 'C2H3IO2'),
       'CH3COOH': (4.76, 'CC(=O)O', 'C2H4O2'), 'Cl2CHCOOH': (1.35, 'ClC(Cl)C(=O)O', 'C2H2Cl2O2'),
       'Cl3CCOOH': (0.66, 'ClC(Cl)(Cl)C(=O)O', 'C2HCl3O2'),
       'A': (2.86, 'CCC(Cl)C(=O)O', 'C4H7ClO2'), 'B': (4.05, 'CC(Cl)CC(=O)O', 'C4H7ClO2'), 'C': (4.52, 'ClCCCC(=O)O', 'C4H7ClO2'),
       'HCOOH': (3.75, 'OC=O', 'CH2O2'), '(CH3)3CCOOH': (5.03, 'CC(C)(C)C(=O)O', 'C5H10O2')}
for k, (p, smi, f) in pKa.items(): assert CalcMolFormula(Chem.MolFromSmiles(smi)) == f, k
order = lambda ks: sorted(ks, key=lambda k: pKa[k][0])
assert order(['ICH2COOH','BrCH2COOH','FCH2COOH','ClCH2COOH']) == ['FCH2COOH','ClCH2COOH','BrCH2COOH','ICH2COOH']
g = ['CH3COOH','ClCH2COOH','Cl2CHCOOH','Cl3CCOOH']; assert order(g)[0] == 'Cl3CCOOH' and order(g)[-1] == 'CH3COOH'
r = 10**(4.76-0.66); print(r); assert sig(r,2) == 1.3e4
assert order(['A','B','C']) == ['A','B','C']
assert order(['CH3COOH','(CH3)3CCOOH','HCOOH']) == ['HCOOH','CH3COOH','(CH3)3CCOOH']

print('全検証OK:', ITEM)
