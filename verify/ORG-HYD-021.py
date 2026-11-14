ITEM = '/home/claude/chem2/out/items/ORG-HYD-021/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
eqs = check_all_equations(); print(len(eqs)); assert len(eqs) == 3
assert check_equation('(CH3)3CCOONa + NaOH -> (CH3)3CH + Na2CO3') is None
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
def salts(alk):
    m = Chem.AddHs(Chem.MolFromSmiles(alk)); out = set()
    for a in m.GetAtoms():
        if a.GetSymbol() != 'H': continue
        rw = Chem.RWMol(m); idx = a.GetIdx(); nb = a.GetNeighbors()[0].GetIdx()
        rw.RemoveAtom(idx)
        c = rw.AddAtom(Chem.Atom(6)); o1 = rw.AddAtom(Chem.Atom(8)); o2 = rw.AddAtom(Chem.Atom(8)); na = rw.AddAtom(Chem.Atom(11))
        nb2 = nb if nb < idx else nb - 1
        rw.AddBond(nb2, c, Chem.BondType.SINGLE); rw.AddBond(c, o1, Chem.BondType.DOUBLE); rw.AddBond(c, o2, Chem.BondType.SINGLE); rw.AddBond(o2, na, Chem.BondType.SINGLE)
        mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(Chem.RemoveHs(mm)))
    return out
p = salts('CCC'); print(p); assert p == {can('CCCC(=O)O[Na]'), can('CC(C)C(=O)O[Na]')}
assert all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H7NaO2' for s in p)
q = salts('CC(C)CC'); print(q)
assert q == {can(s) for s in ['CC(C)CCC(=O)O[Na]', 'CCC(C)CC(=O)O[Na]', 'CC(C)C(C)C(=O)O[Na]', 'CCC(C)(C)C(=O)O[Na]']}

print('全検証OK:', ITEM)
