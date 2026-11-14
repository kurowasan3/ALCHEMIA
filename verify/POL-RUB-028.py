# 検証スクリプト POL-RUB-028  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/POL-RUB-028/item.json'
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


from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
m = Chem.MolFromSmiles('CC1OC(=O)C(C)OC1=O')
iso = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)}
print(iso); assert len(iso) == 3
def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@','!').replace('@','@@').replace('!','@')))
meso = [s for s in iso if mirror(s) == s]; print('meso', meso); assert len(meso) == 1
assert CalcMolFormula(Chem.MolFromSmiles('CC(O)C(=O)O')) == 'C3H6O3'
assert 3*12 + 4*1.0 + 2*16 == 72 and 7.2e4/72 == 1000
print('OK')
