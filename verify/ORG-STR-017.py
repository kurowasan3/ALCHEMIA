ITEM = '/home/claude/chem2/out/items/ORG-STR-017/item.json'
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
assert not check_all_equations()
import itertools
from fractions import Fraction as F
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog('rdApp.*')
def U(m):
    c = h = x = n = 0
    for a in Chem.AddHs(m).GetAtoms():
        s = a.GetSymbol()
        c += s == 'C'; h += s == 'H'; x += s in ('F', 'Cl', 'Br', 'I'); n += s == 'N'
    return F(c) - F(h + x, 2) + F(n, 2) + 1
def rings_pi(m):
    pi = 0
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE: pi += 1
        elif b.GetBondType() == Chem.BondType.TRIPLE: pi += 2
    return m.GetRingInfo().NumRings() + pi
for smi, f, u in [('C=CC#N', 'C3H3N', 3), ('C=CCl', 'C2H3Cl', 1), ('O=[N+]([O-])c1ccccc1', 'C6H5NO2', 5), ('c1ccncc1', 'C5H5N', 4)]:
    m = Chem.MolFromSmiles(smi); Chem.Kekulize(m, clearAromaticFlags=True)
    assert CalcMolFormula(m) == f, (smi, CalcMolFormula(m))
    assert U(m) == u == rings_pi(m), (smi, U(m), rings_pi(m))
    print('OK', f, u)
assert CalcMolFormula(Chem.MolFromSmiles('CNC')) == 'C2H7N'
# (5) C3H6O の鎖式異性体を総当たり
atoms = [6, 6, 6, 8]; pairs = list(itertools.combinations(range(4), 2)); found = set()
for orders in itertools.product(range(3), repeat=len(pairs)):
    rw = Chem.RWMol()
    for z in atoms: rw.AddAtom(Chem.Atom(z))
    for (a, b), o in zip(pairs, orders):
        if o: rw.AddBond(a, b, Chem.BondType.SINGLE if o == 1 else Chem.BondType.DOUBLE)
    m = rw.GetMol()
    try: Chem.SanitizeMol(m)
    except Exception: continue
    if len(Chem.GetMolFrags(m)) != 1 or m.GetRingInfo().NumRings(): continue
    if CalcMolFormula(m) != 'C3H6O': continue
    found.add(Chem.MolToSmiles(m))
print(sorted(found))
noCC = sorted(s for s in found if not Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('C=C')))
print(noCC); assert noCC == sorted([Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in ['CCC=O', 'CC(C)=O']])

print('全検証OK:', ITEM)
