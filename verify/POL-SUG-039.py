# 検証スクリプト（グループ g08）  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/POL-SUG-039/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if should_check(b):
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n_eq += 1
print('反応式のつり合い OK：', n_eq, '式')
# 小問数と答えの数の一致
m_en = re.search(r'\\begin\{enumerate\}(.*)\\end\{enumerate\}', it['problem'], re.S)
nq = m_en.group(1).count('\\item') if m_en else 0
assert '\\begin{itemize}' not in it['problem']
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問数', nq, '答えの数', na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
P = Chem.MolFromSmiles('OCC1OC(Oc2ccc(O)cc2)C(O)C(O)C1O')
Q = Chem.MolFromSmiles('OCC1OC(Oc2ccc(CO)cc2)C(O)C(O)C1O')
Qi = Chem.MolFromSmiles('OCC1OC(OCc2ccc(O)cc2)C(O)C(O)C1O')
print(CalcMolFormula(P), CalcMolFormula(Q), CalcMolFormula(Qi))
assert CalcMolFormula(P)=='C12H16O7' and CalcMolFormula(Q)=='C13H18O7' and CalcMolFormula(Qi)=='C13H18O7'
assert 12*12+16*1+7*16 == 272
ph = Chem.MolFromSmarts('c[OX2H]')
assert P.HasSubstructMatch(ph) and not Q.HasSubstructMatch(ph) and Qi.HasSubstructMatch(ph)
from sympy import Rational as R
C = R(528,10)*12/44; H = R(144,10)*2/18; O = R(272,10)-C-H
r = [C/12, H/1, O/16]; m = min(r); print([x/m for x in r])
assert [x/r[2]*7 for x in r] == [12,16,7]

print('ALL OK')
