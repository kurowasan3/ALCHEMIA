ITEM = '/home/claude/chem2/outh/items/ORG-HYD-022/item.json'
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
from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
RDLogger.DisableLog('rdApp.*')
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
mirror = lambda s: can(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
def clean(rw):
    for a in rw.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    m = rw.GetMol(); Chem.SanitizeMol(m); return Chem.MolToSmiles(m)
def chlorides(alk):
    m = Chem.AddHs(Chem.MolFromSmiles(alk)); out = set()
    for a in m.GetAtoms():
        if a.GetSymbol() == 'H':
            rw = Chem.RWMol(m); rw.GetAtomWithIdx(a.GetIdx()).SetAtomicNum(17); out.add(Chem.MolToSmiles(Chem.RemoveHs(rw)))
    return out
def eliminate(cl):
    m = Chem.MolFromSmiles(cl); out = set()
    clat = [a for a in m.GetAtoms() if a.GetSymbol() == 'Cl'][0]; c = clat.GetNeighbors()[0]
    for nb in c.GetNeighbors():
        if nb.GetSymbol() == 'C' and nb.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            rw.RemoveAtom(clat.GetIdx()); out.add(clean(rw))
    return out
def by_bonds(alk):
    m = Chem.MolFromSmiles(alk); out = set()
    for b in m.GetBonds():
        if b.GetBeginAtom().GetTotalNumHs() and b.GetEndAtom().GetTotalNumHs():
            rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.DOUBLE); out.add(clean(rw))
    return out
def st(s): return sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s), options=StereoEnumerationOptions(unique=True))})
H = {'ヘキサン': 'CCCCCC', '2-メチルペンタン': 'CC(C)CCC', '3-メチルペンタン': 'CCC(C)CC', '2,2-ジメチルブタン': 'CC(C)(C)CC', '2,3-ジメチルブタン': 'CC(C)C(C)C'}
info = {}
for name, s in H.items():
    cls = chlorides(s); al = set().union(*[eliminate(c) for c in cls])
    assert al == by_bonds(s), name
    chiral = [x for x in al if any(mirror(y) != y for y in st(x))]
    ct = [x for x in al if len(st(x)) > 1 and x not in chiral]
    noH = [x for x in al if any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetBeginAtom().GetTotalNumHs() == 0 and b.GetEndAtom().GetTotalNumHs() == 0 for b in Chem.MolFromSmiles(x).GetBonds())]
    dead = [c for c in cls if not eliminate(c)]
    info[name] = dict(ncl=len(cls), nal=len(al), nst=sum(len(st(x)) for x in al), chiral=len(chiral), ct=len(ct), noH=len(noH), dead=dead)
    print(name, info[name])
exp = {'ヘキサン': (3, 5), '2-メチルペンタン': (4, 5), '3-メチルペンタン': (3, 5), '2,2-ジメチルブタン': (1, 1), '2,3-ジメチルブタン': (2, 2)}
assert all((info[k]['nal'], info[k]['nst']) == v for k, v in exp.items())
A = [k for k, v in info.items() if v['ncl'] == 3 and v['nal'] == 1]; assert A == ['2,2-ジメチルブタン']
B = [k for k, v in info.items() if v['nal'] == 3 and v['chiral'] == 1]; assert B == ['3-メチルペンタン']
C = [k for k, v in info.items() if v['nal'] == 4 and v['ct'] == 1]; assert C == ['2-メチルペンタン']
D = [k for k, v in info.items() if v['nal'] == 2 and v['noH'] == 1]; assert D == ['2,3-ジメチルブタン']
assert info['2,2-ジメチルブタン']['dead'] == [can('CCC(C)(C)CCl')]
assert info['ヘキサン']['chiral'] == 0 and info['ヘキサン']['ct'] == 2 and info['3-メチルペンタン']['ct'] == 1

print('全検証OK:', ITEM)
