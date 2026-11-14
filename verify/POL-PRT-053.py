# 検証スクリプト POL-PRT-053  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/POL-PRT-053/item.json'
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
AW = {'C':12.0,'H':1.0,'N':14.0,'O':16.0}
def mw(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); return sum(AW[a.GetSymbol()] for a in m.GetAtoms())
for s, v in (('NCC(=O)O',75),('CC(N)C(=O)O',89),('OCC(N)C(=O)O',105),('NC(Cc1ccccc1)C(=O)O',165)):
    assert mw(s) == v
cyc = 'O=C1CNC(=O)CNC(=O)C(C)NC(=O)C(C)NC(=O)C(CO)NC(=O)C(Cc2ccccc2)N1'
lin = 'NCC(=O)NCC(=O)NC(C)C(=O)NC(C)C(=O)NC(CO)C(=O)NC(Cc1ccccc1)C(=O)O'
print(mw(cyc), mw(lin)); assert mw(cyc) == 490 and mw(lin) == 508
w = 0.98/490*6*18; print(w); assert sig(w,2) == 0.22
pat = Chem.MolFromSmarts('[NX3;H2]')
assert not Chem.MolFromSmiles(cyc).HasSubstructMatch(pat) and Chem.MolFromSmiles(lin).HasSubstructMatch(pat)
print('OK')
