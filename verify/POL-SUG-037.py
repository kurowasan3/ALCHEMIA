# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/POL-SUG-037/item.json'
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

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from fractions import Fraction
AM = dict(H=1.0, C=12.0, N=14.0, O=16.0, Na=23.0, K=39.0, Cu=63.5, I=127.0)
from chemlib import parse_formula
def M(f):
    c, _ = parse_formula(f)
    return sum(AM[e] * float(n) for e, n in c.items())

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem import rdCIPLabeler
# D-グルコース（鎖状）を立体異性体の列挙から CIP (2R,3S,4R,5R) で特定する
base = Chem.MolFromSmiles('O=CC(O)C(O)C(O)C(O)CO')
def cip(m):
    rdCIPLabeler.AssignCIPLabels(m)
    lab = {}
    for a in m.GetAtoms():
        if a.HasProp('_CIPCode'): lab[a.GetIdx()] = a.GetProp('_CIPCode')
    return lab
# 原子番号: 0 O,1 C1,2 C2,3 O,4 C3,5 O,6 C4,7 O,8 C5,9 O,10 C6,11 O
CI = {2: 2, 3: 4, 4: 6, 5: 8}
glc = None
for iso in EnumerateStereoisomers(base):
    l = cip(Chem.Mol(iso))
    if [l[CI[k]] for k in (2,3,4,5)] == ['R','S','R','R']:
        glc = iso
assert glc is not None
GLC = 'RLRR'
ALD = dict(アロース='RRRR', アルトロース='LRRR', グルコース='RLRR', マンノース='LLRR', グロース='RRLR', イドース='LRLR', ガラクトース='RLLR', タロース='LLLR')
KET = dict(プシコース='RRR', フルクトース='LRR', ソルボース='RLR', タガトース='LLR')
def aldose(pat):
    m = Chem.RWMol(glc)
    for k, c in zip((2,3,4,5), pat):
        if c != GLC[k-2]:
            a = m.GetAtomWithIdx(CI[k]); a.InvertChirality()
    return m
def can(m): return Chem.MolToSmiles(m)
def mirror(m):
    m = Chem.RWMol(m)
    for a in m.GetAtoms():
        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED: a.InvertChirality()
    return m
def alditol(m):
    m = Chem.RWMol(m); m.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
    m.UpdatePropertyCache(False); Chem.SanitizeMol(m); return m
def ketose(m):  # C1 を CH2OH，C2 を C=O
    m = Chem.RWMol(m); m.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
    m.GetAtomWithIdx(2).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    m.GetBondBetweenAtoms(2, 3).SetBondType(Chem.BondType.DOUBLE)
    m.UpdatePropertyCache(False); Chem.SanitizeMol(m); return m
def enediol12(m):
    m = Chem.RWMol(m); m.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
    m.GetAtomWithIdx(2).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    m.GetBondBetweenAtoms(1, 2).SetBondType(Chem.BondType.DOUBLE)
    m.UpdatePropertyCache(False); Chem.SanitizeMol(m); return m
A = {k: aldose(v) for k, v in ALD.items()}
# 鏡像の関係・D 系列の確認：CIP で D-グルコース
# ケトースの名前：C3〜C5 のパターンで対応づけ
K = {}
for k, v in ALD.items():
    kn = [n for n, p in KET.items() if p == v[1:]][0]
    kk = can(ketose(A[k]))
    K.setdefault(kn, kk); assert K[kn] == kk
# (1) グルコースと同じ 1,2-エンジオール
e = can(enediol12(A['グルコース']))
same = [k for k in A if can(enediol12(A[k])) == e]
assert sorted(same) == sorted(['グルコース', 'マンノース'])
assert can(ketose(A['グルコース'])) == K['フルクトース']
# (2) 還元生成物がメソ体（光学不活性）のアルドース
meso = [k for k in A if can(alditol(A[k])) == can(mirror(alditol(A[k])))]
assert sorted(meso) == sorted(['アロース', 'ガラクトース']), meso
cands = []
for p in meso:
    for q in A:
        if q != p and can(enediol12(A[q])) == can(enediol12(A[p])):
            assert can(alditol(A[q])) != can(mirror(alditol(A[q])))   # Q は光学活性
            kn = [n for n, kk in K.items() if kk == can(ketose(A[p]))][0]
            cands.append((p, q, kn))
print('(2) 候補', cands)
assert sorted(cands) == sorted([('アロース','アルトロース','プシコース'), ('ガラクトース','タロース','タガトース')])
# (3) C3 の配置を反転した 2-ケトース
flip3 = lambda p: ('L' if p[0]=='R' else 'R') + p[1:]
ok = [c for c in cands if flip3(KET[c[2]]) == KET['フルクトース']]
assert ok == [('アロース','アルトロース','プシコース')]
# (4) 到達できる糖（C2, C3 だけ変わる）
reach_a = {k for k, v in ALD.items() if v[2:] == ALD['グルコース'][2:]}
reach_k = {k for k, v in KET.items() if v[1:] == KET['フルクトース'][1:]}
assert reach_a == {'アロース','アルトロース','グルコース','マンノース'} and reach_k == {'プシコース','フルクトース'}
# L-グルコースは D-グルコースの鏡像（C4, C5 も異なる）
assert can(mirror(A['グルコース'])) != can(A['グルコース'])
print('POL-SUG-037 構造OK')
