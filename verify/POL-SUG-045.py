# 検証スクリプト POL-SUG-045  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/POL-SUG-045/item.json'
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
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


from fractions import Fraction as Fr
C = 0.792*12/44; H = 0.288*2.0/18; O = 0.504 - C - H
r = [C/12, H/1.0, O/16]; m = min(r); print([x/m for x in r])
assert [round(x/(r[0]/9)) for x in r] == [9, 16, 8]
M = 18*12 + 32*1.0 + 16*16; assert M == 504
assert abs(0.504/M*18*44 - 0.792) < 1e-9 and abs(0.504/M*16*18 - 0.288) < 1e-9
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
glc = 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'
suc = 'OC[C@H]1O[C@H](O[C@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O'
f = [CalcMolFormula(Chem.MolFromSmiles(s)) for s in (glc, suc)]
print(f); assert f == ['C6H12O6', 'C12H22O11']
ratios = {Fr(6,12), Fr(12,22), Fr(18,32)}; assert len(ratios) == 3
print('OK')
