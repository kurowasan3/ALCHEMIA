# 検証スクリプト POL-RUB-019  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/POL-RUB-019/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
from fractions import Fraction as F
for n in (1, 2, 5):
    assert check_equation(f'C{3*n}H{4*n}O{2*n} + {3*n}O2 -> {3*n}CO2 + {2*n}H2O') is None
    assert check_equation(f'{n}C3H6O3 -> C{3*n}H{4*n}O{2*n} + {n}H2O') is None
u = F(144,10)/72
assert u == F(2,10)
co2 = 3*u; v = co2*F(224,10)
assert sig(float(v), 3) == 13.4 and 3*u == F(6,10)
glc = co2/6*180; pla = co2/3*72
assert glc == 18 and pla == F(144,10)
pe = F(144,10)/14*F(224,10)
assert sig(float(pe), 3) == 23.0, float(pe)
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert CalcMolFormula(Chem.MolFromSmiles('CC(O)C(=O)O')) == 'C3H6O3'
from collections import Counter
cnt = Counter(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles('*OC(C)C(=O)*')).GetAtoms() if a.GetSymbol()!='*')
assert cnt == Counter({'C':3,'H':4,'O':2})
print('炭素の循環 OK', float(v), float(pe))

print('ALL OK')
