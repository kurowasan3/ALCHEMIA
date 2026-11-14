# ORG-STR-012 検証：C6H12O（鎖式・エノール除く）の異性体を列挙し，条件で A〜D が一意に決まるか，立体異性体数を確認
import json, re, sys
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def canon(m):
    try:
        mm = m.GetMol() if hasattr(m, 'GetMol') else m
        Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
    except Exception:
        return None
alk = ['CCCCCC', 'CC(C)CCC', 'CCC(C)CC', 'CC(C)C(C)C', 'CC(C)(C)CC']
sat = set()
for s in alk:
    m = Chem.MolFromSmiles(s)
    for a in m.GetAtoms():
        rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(8)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE); sat.add(canon(rw))
    for b in m.GetBonds():
        rw = Chem.RWMol(m); i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); rw.RemoveBond(i, j)
        k = rw.AddAtom(Chem.Atom(8)); rw.AddBond(i, k, Chem.BondType.SINGLE); rw.AddBond(k, j, Chem.BondType.SINGLE); sat.add(canon(rw))
sat.discard(None)
assert len(sat) == 32   # C6H14O：アルコール17＋エーテル15
uns = set()
for s in sat:
    m = Chem.MolFromSmiles(s)
    for b in m.GetBonds():
        if b.GetBeginAtom().GetTotalNumHs() and b.GetEndAtom().GetTotalNumHs():
            rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.DOUBLE); uns.add(canon(rw))
uns.discard(None)
enol = Chem.MolFromSmarts('C=C[OX2H]')
L = sorted(s for s in uns if not Chem.MolFromSmiles(s).HasSubstructMatch(enol))
for s in L:
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C6H12O'
P = Chem.MolFromSmarts
OH, CC, ALD, KET = P('[CX4][OX2H]'), P('C=C'), P('[CX3H1]=O'), P('[#6][CX3](=O)[#6]')
IOD = [P('[CH3]C(=O)'), P('[CH3][CX4H1][OX2H]')]
def nchir(m):
    si = Chem.FindPotentialStereo(Chem.AddHs(m))
    return sum(1 for x in si if x.type == Chem.StereoType.Atom_Tetrahedral)
def nct(m):
    return sum(1 for x in Chem.FindPotentialStereo(m) if x.type == Chem.StereoType.Bond_Double)
def hyd(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s))
    for b in m.GetBonds(): b.SetBondType(Chem.BondType.SINGLE)
    return canon(m)
def ozon(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s))
    b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE and b.GetBeginAtom().GetSymbol() == b.GetEndAtom().GetSymbol() == 'C'][0]
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); m.RemoveBond(i, j)
    for k in (i, j):
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(k, o, Chem.BondType.DOUBLE)
    mm = m.GetMol(); Chem.SanitizeMol(mm)
    return [Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mm, asMols=True)]
def ohdeg(m):
    c = m.GetSubstructMatch(OH)[0]; return sum(1 for n in m.GetAtomWithIdx(c).GetNeighbors() if n.GetSymbol() == 'C')
M = {s: Chem.MolFromSmiles(s) for s in L}
alc = [s for s in L if M[s].HasSubstructMatch(OH)]
assert all(M[s].HasSubstructMatch(CC) for s in alc)            # C6H12O の鎖式アルコールはすべて C=C をもつ
noNa_noBr = [s for s in L if not M[s].HasSubstructMatch(OH) and not M[s].HasSubstructMatch(CC)]
assert all(M[s].HasSubstructMatch(ALD) or M[s].HasSubstructMatch(KET) for s in noNa_noBr)
D = [s for s in noNa_noBr if not M[s].HasSubstructMatch(ALD) and nchir(M[s]) >= 1 and any(M[s].HasSubstructMatch(p) for p in IOD)]
Apre = [s for s in alc if nct(M[s]) >= 1 and nchir(M[s]) == 0 and nchir(Chem.MolFromSmiles(hyd(s))) >= 1]
A = [s for s in Apre if any(Chem.MolFromSmiles(p).HasSubstructMatch(KET) and not Chem.MolFromSmiles(p).HasSubstructMatch(P('[OX2H]')) for p in ozon(s))]
B = [s for s in alc if nct(M[s]) == 0 and nchir(M[s]) >= 1 and nchir(Chem.MolFromSmiles(hyd(s))) == 0 and ohdeg(M[s]) == 3]
Bpre = [s for s in alc if nct(M[s]) == 0 and nchir(M[s]) >= 1 and nchir(Chem.MolFromSmiles(hyd(s))) == 0]
Cpre = [s for s in alc if nct(M[s]) >= 1 and nchir(M[s]) >= 1 and any(M[s].HasSubstructMatch(p) for p in IOD)]
C = [s for s in Cpre if nchir(Chem.MolFromSmiles(hyd(s))) == 2]
print('isomers', len(L)); print('Apre', [(s, ozon(s)) for s in Apre]); print('A', A, 'Bpre', Bpre, 'B', B, 'Cpre', Cpre, 'C', C, 'D', D)
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
assert Apre == sorted(map(can, ['CC=C(C)CCO', 'CCC(C)=CCO', 'CCC=C(C)CO']))
assert A == [can('CCC(C)=CCO')] and B == [can('C=CC(C)(O)CC')] and C == [can('CC=C(C)C(C)O')] and D == [can('CCC(C)C(C)=O')]
assert sorted(Bpre) == sorted(map(can, ['C=CC(C)(O)CC', 'C=CC(CC)CO']))
assert sorted(Cpre) == sorted(map(can, ['CC=C(C)C(C)O', 'CC=CCC(C)O', 'CCC=CC(C)O']))
G = hyd(C[0]); Dred = hyd(D[0]); assert G == Dred == can('CCC(C)C(C)O')
nst = lambda s: len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s))})
assert nst(C[0]) == 4 and nst(G) == 4
assert nchir(Chem.MolFromSmiles(hyd(A[0]))) == 1 and nchir(Chem.MolFromSmiles(hyd(B[0]))) == 0
# (4) C の1立体異性体 → G：C2 の配置は保持，C3 は両方 → 2 種，互いに鏡像でない
c1 = 'C/C=C(C)/[C@@H](C)O'
prods = {Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in ['CC[C@H](C)[C@@H](C)O', 'CC[C@@H](C)[C@@H](C)O']}
mirror = lambda s: s.replace('@@', '!').replace('@', '@@').replace('!', '@')
p = sorted(prods); assert len(p) == 2 and Chem.MolToSmiles(Chem.MolFromSmiles(mirror(p[0]))) != p[1]
# 反応式の確認（item.json の \ce{} で矢印1つの完全な式）
it = json.load(open('/home/claude/chem/out2/items/ORG-STR-012/item.json'))
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') == 1 and 'boxed' not in eq and 'R' not in eq:
            r = check_equation(eq); assert r is None, (eq, r); print('eq ok', eq)
print('OK')
