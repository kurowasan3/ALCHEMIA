# 検証スクリプト POL-RUB-020  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/POL-RUB-020/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
from sympy import Rational as R, symbols, solve
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
check_structure(); check_all_equations()
# 原子量 H=1, C=12, O=16, Br=80 での式量
def fw(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); w = {'H':1,'C':12,'O':16,'Br':80}
    return sum(w[a.GetSymbol()] for a in m.GetAtoms())
unit = fw('CC(O)')-2        # -CH2-CH(OH)-  (両端の H を除く) = C2H4O
cin = fw('OC(=O)/C=C/c1ccccc1'); water = fw('O')
est = unit + cin - water
print('単位', unit, 'ケイ皮酸', cin, 'エステル単位', est)
assert (unit, cin, est) == (44, 148, 174)
n = symbols('n'); nn = solve(R(220,10) + 130*n - R(610,10), n)[0]
assert nn == R(3,10); frac = nn/(R(22)/44); assert frac == R(3,5)
br_before = R(3,100)*160; assert br_before == R(48,10)
rem = R(288,100)/160; reacted = R(3,100)-rem
print('残 C=C', rem, '反応', reacted, '割合', reacted/R(3,100), '四員環', reacted/2)
assert rem == R(18,1000) and reacted/R(3,100) == R(2,5) and reacted/2 == R(6,1000)
# 二量化で原子の出入りなし
mono = Chem.MolFromSmiles('COC(=O)/C=C/c1ccccc1'); dim = Chem.MolFromSmiles('COC(=O)C1C(c2ccccc2)C(C(=O)OC)C1c1ccccc1')
print(CalcMolFormula(mono), CalcMolFormula(dim)); assert abs(2*MolWt(mono)-MolWt(dim)) < 1e-6
print('POL-RUB-020 OK')
