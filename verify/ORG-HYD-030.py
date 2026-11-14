# 検証スクリプト ORG-HYD-030（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem5/out/items/ORG-HYD-030/item.json')
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
def can(m): return Chem.MolToSmiles(m)
def alkanes(n):
    cur={'C'}
    for _ in range(n-1):
        nxt=set()
        for s in cur:
            m=Chem.MolFromSmiles(s)
            for a in m.GetAtoms():
                if a.GetDegree()<4:
                    rw=Chem.RWMol(m); i=rw.AddAtom(Chem.Atom(6)); rw.AddBond(a.GetIdx(),i,Chem.BondType.SINGLE)
                    nxt.add(can(rw.GetMol()))
        cur=nxt
    return cur
def alkenes(n):
    out=set()
    for s in alkanes(n):
        m=Chem.MolFromSmiles(s)
        for b in m.GetBonds():
            a1,a2=b.GetBeginAtom(),b.GetEndAtom()
            if a1.GetTotalNumHs()>0 and a2.GetTotalNumHs()>0:
                rw=Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
                mm=rw.GetMol(); Chem.SanitizeMol(mm); out.add(can(mm))
    return out
def cycloalkanes(n):
    out=set()
    for s in alkanes(n):
        m=Chem.MolFromSmiles(s)
        N=m.GetNumAtoms()
        for i in range(N):
            for j in range(i+1,N):
                if m.GetBondBetweenAtoms(i,j): continue
                if m.GetAtomWithIdx(i).GetDegree()<4 and m.GetAtomWithIdx(j).GetDegree()<4:
                    rw=Chem.RWMol(m); rw.AddBond(i,j,Chem.BondType.SINGLE); mm=rw.GetMol(); Chem.SanitizeMol(mm)
                    out.add(can(mm))
    return out
def stereo(s):
    m=Chem.MolFromSmiles(s)
    iso=[Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m,options=StereoEnumerationOptions(unique=True))]
    iso=sorted(set(iso))
    def mirror(x):
        return Chem.MolToSmiles(Chem.MolFromSmiles(x.replace('@@','!').replace('@','@@').replace('!','@')))
    chiral=[x for x in iso if mirror(x)!=x]
    return iso, chiral

def add(s, X, Y):
    m = Chem.MolFromSmiles(s); prods = set()
    for b in m.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE: continue
        for (p, q) in [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), (b.GetEndAtomIdx(), b.GetBeginAtomIdx())]:
            rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.SINGLE)
            i = rw.AddAtom(Chem.Atom(X)); rw.AddBond(p, i, Chem.BondType.SINGLE)
            if Y:
                j = rw.AddAtom(Chem.Atom(Y)); rw.AddBond(q, j, Chem.BondType.SINGLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); prods.add(can(mm))
    return prods
nst = lambda s: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
first = {}
for n in range(2, 7):
    for s in sorted(alkenes(n)):
        if stereo(s)[1]: first.setdefault(1, (n, set())); first[1][1].add(s) if first[1][0] == n else None
        if any(nst(p) >= 1 for p in add(s, 35, 35)): first.setdefault(2, (n, set())); first[2][1].add(s) if first[2][0] == n else None
        if len(add(s, 8, 0)) >= 2: first.setdefault(3, (n, set())); first[3][1].add(s) if first[3][0] == n else None
        if any(nst(p) >= 2 for p in add(s, 35, 35)): first.setdefault(4, (n, set())); first[4][1].add(s) if first[4][0] == n else None
for n in range(3, 7):
    for s in sorted(cycloalkanes(n)):
        if stereo(s)[1]: first.setdefault(5, (n, set())); first[5][1].add(s) if first[5][0] == n else None
print(first)
assert first[1] == (6, {'C=CC(C)CC'})
assert first[2] == (3, {'C=CC'}) and first[3] == (3, {'C=CC'})
assert first[4] == (4, {'CC=CC'})
assert first[5] == (5, {'CC1CC1C'})
c4one = [s for s in alkenes(4) if len(add(s, 8, 0)) == 1]
assert c4one == ['CC=CC']

print('ALL OK')
