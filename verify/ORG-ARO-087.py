# 検証スクリプト（東工大 g12）  python3 -I この.py
import json, re, itertools
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
RDLogger.DisableLog('rdApp.*')
ITEM = '/home/claude/chemx/tokodai/outh/items/ORG-ARO-087/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
assert it['problem'].count('\\item') == len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])) == 5
# 分子量
assert 226 + 60 - 18 == 268 and 16*12 + 18 + 16 == 226
# 候補の列挙：C4 の飽和骨格（直鎖・枝分かれ）に Ph 2 個と OH 1 個（全員炭素に結合）
skels = {'n': [(0,1),(1,2),(2,3)], 'iso': [(0,1),(1,2),(1,3)]}
def build(bonds, phs, oh):
    m = Chem.RWMol()
    for _ in range(4): m.AddAtom(Chem.Atom(6))
    for a,b in bonds: m.AddBond(a,b,Chem.BondType.SINGLE)
    o = m.AddAtom(Chem.Atom(8)); m.AddBond(oh,o,Chem.BondType.SINGLE)
    smi = Chem.MolToSmiles(m.GetMol())
    s = Chem.MolFromSmiles(smi)
    return s
cands = set()
for bonds in skels.values():
    for phs in itertools.combinations_with_replacement(range(4), 2):
        for oh in range(4):
            m = Chem.RWMol()
            for _ in range(4): m.AddAtom(Chem.Atom(6))
            for a,b in bonds: m.AddBond(a,b,Chem.BondType.SINGLE)
            o = m.AddAtom(Chem.Atom(8)); m.AddBond(oh,o,Chem.BondType.SINGLE)
            for p in phs:
                ring = [m.AddAtom(Chem.Atom(6)) for _ in range(6)]
                for i in range(6): m.AddBond(ring[i], ring[(i+1)%6], Chem.BondType.AROMATIC)
                for r in ring: m.GetAtomWithIdx(r).SetIsAromatic(True)
                m.AddBond(p, ring[0], Chem.BondType.SINGLE)
            try:
                mol = m.GetMol(); Chem.SanitizeMol(mol); cands.add(Chem.MolToSmiles(mol))
            except Exception: pass
for s in cands: assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C16H18O'
print('ア・イを満たす構造異性体', len(cands))
def cls(s):
    m = Chem.MolFromSmiles(s); o = [a for a in m.GetAtoms() if a.GetSymbol()=='O'][0]
    return sum(1 for n in o.GetNeighbors()[0].GetNeighbors() if n.GetSymbol()=='C')
def nst(s): return len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s))})
def nchiral(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
def nEZ(s): return sum(1 for x in Chem.FindPotentialStereo(Chem.MolFromSmiles(s)) if x.type == Chem.StereoType.Bond_Double)
dehyd = AllChem.ReactionFromSmarts('[C;!H0;X4:1]-[C;X4:2]-[OX2H]>>[C:1]=[C:2]')
h2 = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]-[C:2]')
def run(r, s):
    out = set()
    for p in r.RunReactants((Chem.MolFromSmiles(s),)):
        q = p[0]
        try: Chem.SanitizeMol(q); out.add(Chem.MolToSmiles(q))
        except Exception: pass
    return sorted(out)
U = lambda s: cls(s) == 2                       # ウ 第二級
E = lambda s: nst(s) == 4                       # エ
def O(s):                                       # オ
    al = run(dehyd, s); return len(al) == 1 and nEZ(al[0]) == 1
def K(s):                                       # カ（B が1種類のときのみ意味をもつ）
    al = run(dehyd, s); return len(al) >= 1 and all(nchiral(run(h2, a)[0]) == 0 for a in al)
sec = sorted(s for s in cands if U(s)); print('ア〜ウ', len(sec)); assert len(sec) == 11
full = [s for s in cands if U(s) and E(s) and O(s) and K(s)]; print('ア〜カ', full)
A = Chem.MolToSmiles(Chem.MolFromSmiles('CC(Cc1ccccc1)C(O)c1ccccc1'))
assert full == [A]
noK = sorted(s for s in cands if U(s) and E(s) and O(s)); print('カなし', noK)
assert set(noK) == {Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in ['CC(Cc1ccccc1)C(O)c1ccccc1','CC(CC(O)c1ccccc1)c1ccccc1','CCC(c1ccccc1)C(O)c1ccccc1']}
# 各条件の必要性
r0=[s for s in cands if E(s) and O(s) and K(s)]; print('ウなし', r0); assert r0==[A]
r = sorted(s for s in cands if U(s) and O(s) and K(s)); print('エなし', r)
assert r == sorted([A, Chem.MolToSmiles(Chem.MolFromSmiles('OC(CCCc1ccccc1)c1ccccc1'))])
# B, C
B = run(dehyd, A); assert len(B) == 1 and nEZ(B[0]) == 1
C = run(h2, B[0])[0]; assert nchiral(C) == 0
assert B[0] == Chem.MolToSmiles(Chem.MolFromSmiles('CC(=Cc1ccccc1)Cc1ccccc1'))
assert C == Chem.MolToSmiles(Chem.MolFromSmiles('CC(Cc1ccccc1)Cc1ccccc1'))
assert nst(A) == 4 and nchiral(A) == 2
# 解説の表（11 種）の判定を再現：立体異性体の数とアルケンの数
table = {'CC(C)(c1ccccc1)C(O)c1ccccc1':(2,0),'CC(CC(O)c1ccccc1)c1ccccc1':(4,1),'CC(Cc1ccccc1)C(O)c1ccccc1':(4,1),
 'CC(O)C(C)(c1ccccc1)c1ccccc1':(2,1),'CC(O)C(Cc1ccccc1)c1ccccc1':(4,2),'CC(O)CC(c1ccccc1)c1ccccc1':(2,2),
 'CC(c1ccccc1)C(O)Cc1ccccc1':(4,2),'CCC(O)C(c1ccccc1)c1ccccc1':(2,2),'CCC(c1ccccc1)C(O)c1ccccc1':(4,1),
 'OC(CCCc1ccccc1)c1ccccc1':(2,1),'OC(CCc1ccccc1)Cc1ccccc1':(2,2)}
for k,(n,a) in table.items():
    s = Chem.MolToSmiles(Chem.MolFromSmiles(k)); assert s in sec, k
    assert nst(s) == n and len(run(dehyd, s)) == a, (k, nst(s), run(dehyd, s))
print('ALL OK')
