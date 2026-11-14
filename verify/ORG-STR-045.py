# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/ORG-STR-045/item.json'
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
assert np_ == na, ('小問の数', np_, na)

# ---- 計算
from sympy import Rational as R, symbols, solve
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import rdMolDescriptors as D
alk = ['C=CCC', 'CC=CC', 'C=C(C)C']
S = set()
for a in alk:
    m = Chem.MolFromSmiles(a)
    # OH on sp3 carbon
    for at in m.GetAtoms():
        if at.GetTotalNumHs() > 0 and not any(b.GetBondTypeAsDouble() == 2 for b in at.GetBonds()):
            rw = Chem.RWMol(m); o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(at.GetIdx(), o, Chem.BondType.SINGLE)
            q = rw.GetMol(); Chem.SanitizeMol(q); S.add(Chem.MolToSmiles(q))
    # O inserted into C-C single bond
    for b in m.GetBonds():
        if b.GetBondTypeAsDouble() == 1:
            rw = Chem.RWMol(m); i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            rw.RemoveBond(i, j); o = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(i, o, Chem.BondType.SINGLE); rw.AddBond(j, o, Chem.BondType.SINGLE)
            q = rw.GetMol(); Chem.SanitizeMol(q); S.add(Chem.MolToSmiles(q))
S = sorted(S)
def nstereo(m):
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True)
    return len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m, options=opts)})
def br2(smi):
    m = Chem.MolFromSmiles(smi)
    for b in m.GetBonds():
        if b.GetBondTypeAsDouble() == 2:
            rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.SINGLE)
            for k in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()):
                x = rw.AddAtom(Chem.Atom(35)); rw.AddBond(k, x, Chem.BondType.SINGLE)
            q = rw.GetMol(); Chem.SanitizeMol(q); return q
def nchiral(m):
    return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
    tot = 0
    for s in S:
        m = Chem.MolFromSmiles(s)
        assert D.CalcMolFormula(m) == 'C4H8O'
        ns = nstereo(m); tot += ns
        a = br2(s)

# ---- (1)(2)
tot = sum(nstereo(Chem.MolFromSmiles(s)) for s in S)
print('構造異性体', len(S), '立体込み', tot)
assert len(S) == 8 and tot == 11
two = sorted(s for s in S if nchiral(br2(s)) == 2)
print('付加体の不斉炭素2個:', two)
assert set(two) == {Chem.CanonSmiles(x) for x in ['CC=CCO', 'C=CC(C)O', 'CC=COC']}
# ---- (3) P：エーテルで二重結合炭素に O → 加水分解（エノール→カルボニル）
iodo = [Chem.MolFromSmarts('[CH3][CX3](=O)[#6,#1]'), Chem.MolFromSmarts('[CH3][CX3H1]=O'),
        Chem.MolFromSmarts('[CH3][CX4H1]([OH1])'), Chem.MolFromSmarts('[CH3][CH2][OH1]')]
def iodoform(smi):
    m = Chem.MolFromSmiles(smi)
    return any(m.HasSubstructMatch(p) for p in iodo)
assert iodoform('CCO') and iodoform('CC=O') and iodoform('CC(C)=O') and not iodoform('CO') and not iodoform('CCC=O')
def hydrolyze(smi):
    m = Chem.MolFromSmiles(smi)
    if m.HasSubstructMatch(Chem.MolFromSmarts('[OX2H1]')): return None
    patt = Chem.MolFromSmarts('[C:1]=[C:2][O:3][C:4]')
    rx = AllChem.ReactionFromSmarts('[C:1]=[C:2][O:3][C:4]>>[C:1][C:2]=[O:3].[O][C:4]')
    if not m.HasSubstructMatch(patt): return 'no'
    p = rx.RunReactants((m,))[0]
    out = []
    for q in p:
        Chem.SanitizeMol(q); out.append(Chem.MolToSmiles(q))
    return sorted(out)
from rdkit.Chem import AllChem
Pc = []
for s in S:
    h = hydrolyze(s)
    print(s, '加水分解:', h)
    if isinstance(h, list) and all(iodoform(x) for x in h): Pc.append(s)
assert Pc == [Chem.CanonSmiles('C=COCC')], Pc
Qc = [s for s in S if Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[OX2H1]')) and nchiral(Chem.MolFromSmiles(s)) == 1]
assert Qc == [Chem.CanonSmiles('C=CC(C)O')], Qc
# ---- (4)(5) 一方の鏡像異性体 Q からの付加体
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
def prods(qsmi):
    a = br2(qsmi)
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True)
    return {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(a, options=opts)}
pR = prods('C=C[C@@H](C)O'); pS = prods('C=C[C@H](C)O')
print('R 体から:', pR, ' S 体から:', pS)
assert len(pR) == 2 and len(pS) == 2 and not (pR & pS)
def mirror(smi):
    return Chem.CanonSmiles(smi.replace('@@', '!').replace('@', '@@').replace('!', '@'))
assert all(mirror(x) not in pR for x in pR)          # R 体からの生成物どうしは鏡像でない
assert {mirror(x) for x in pR} == pS                  # ラセミ体では鏡像の組が 2 組
print('ラセミ体からの立体異性体', len(pR | pS), '留分', len(pR))
# 反応式（加水分解）
for e in ['CH2=CHOCH2CH3 + H2O -> CH3CHO + CH3CH2OH', 'CH3CH=CHOCH3 + H2O -> CH3CH2CHO + CH3OH', 'CH2=C(CH3)OCH3 + H2O -> CH3COCH3 + CH3OH']:
    assert check_equation(e) is None
print('OK')
