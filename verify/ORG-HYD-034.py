# 検証スクリプト ORG-HYD-034（グループ g01）  python3 -I この.py
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

it = check_item('/home/claude/chem5/outh/items/ORG-HYD-034/item.json')
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

from itertools import combinations
def cleave(s):
    m=Chem.MolFromSmiles(s)
    rw=Chem.RWMol(m)
    dbs=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in m.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE]
    for a,b in dbs:
        rw.RemoveBond(a,b)
        for x in (a,b):
            o=rw.AddAtom(Chem.Atom(8)); rw.AddBond(x,o,Chem.BondType.DOUBLE)
    mm=rw.GetMol(); Chem.SanitizeMol(mm)
    frags=[]
    for f in Chem.GetMolFrags(mm,asMols=True):
        # oxidise aldehyde C-H to COOH
        changed=True
        while changed:
            changed=False
            for a in f.GetAtoms():
                if a.GetSymbol()=='C' and a.GetTotalNumHs()>0 and any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetSymbol()=='O' for b in a.GetBonds()):
                    r=Chem.RWMol(f); o=r.AddAtom(Chem.Atom(8)); r.AddBond(a.GetIdx(),o,Chem.BondType.SINGLE)
                    f=r.GetMol(); Chem.SanitizeMol(f); changed=True; break
        s2=Chem.MolToSmiles(f)
        if s2 in ('O=CO','OC(=O)O','O=C(O)O','C=O'): s2='O=C=O'
        frags.append(s2)
    return tuple(sorted(frags))
def place_db(s,k):
    m=Chem.MolFromSmiles(s); out=set()
    for bs in combinations([b.GetIdx() for b in m.GetBonds()],k):
        rw=Chem.RWMol(m)
        try:
            for bi in bs: rw.GetBondWithIdx(bi).SetBondType(Chem.BondType.DOUBLE)
            mm=rw.GetMol(); Chem.SanitizeMol(mm)
        except Exception: continue
        # no cumulated
        if any(sum(1 for b in a.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE)>1 for a in mm.GetAtoms()): continue
        out.add(can(mm))
    return out

from rdkit.Chem.rdMolDescriptors import CalcMolFormula
C = 99.0 * 12 / 44; Hm = 36.0 * 2.0 / 18
assert abs(C + Hm - 31.0) < 1e-9 and abs((C / 12) / (Hm / 1.0) - 9 / 16) < 1e-9 and 9 * 12 + 16 * 1.0 == 124
target = tuple(sorted(['CC(C)=O', 'O=C=O', 'CC(=O)CCC(=O)O']))
dienes = set()
for s in alkanes(9): dienes |= place_db(s, 2)
hits = sorted(s for s in dienes if cleave(s) == target)
print('A 候補', hits); assert hits == sorted(['C=CCCC(C)=C(C)C', 'C=C(C)CCC=C(C)C'])
achiral = [h for h in hits if not Chem.FindMolChiralCenters(Chem.MolFromSmiles(h.replace('=', '')), includeUnassigned=True, useLegacyImplementation=False)]
assert achiral == ['C=C(C)CCC=C(C)C']
ring = set()
for s in cycloalkanes(9): ring |= place_db(s, 1)
cands = []
for s in ring:
    p = cleave(s)
    if len(p) != 1: continue
    pm = Chem.MolFromSmiles(p[0])
    if CalcMolFormula(pm) != 'C9H16O3': continue
    if len(pm.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[OH]'))) != 1 or not pm.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=O)([#6])[#6]')): continue
    if any(sum(1 for nb in a.GetNeighbors() if nb.GetSymbol() == 'C') > 2 for a in pm.GetAtoms() if a.GetSymbol() == 'C'): continue
    cands.append((s, p[0], pm.HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)[#6]'))))
print(len(cands), cands)
assert len(cands) == 6 and [c[0] for c in cands if c[2]] == ['CC1=CCCCCCC1']

print('ALL OK')
