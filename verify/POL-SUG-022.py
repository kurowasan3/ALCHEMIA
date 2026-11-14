import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def load(item_id):
    return json.load(open(f'/home/claude/chem/out2/items/{item_id}/item.json', encoding='utf-8'))

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

def check_all_equations(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return len(eqs)

def check_subitems(it):
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    if n_items:
        assert [int(x) for x in nums] == list(range(1, n_items + 1)), (n_items, nums)
    print('小問数OK:', n_items)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def in_text(val_str, it, key='answer'):
    assert val_str in it[key], (val_str, key)
# ---- POL-SUG-022 過ヨウ素酸酸化のシミュレーション（RDKit のグラフ操作）
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
def periodate(smi):
    """隣接ジオール（C-OH が隣り合う C-C 結合）をすべて切断．切れた結合 1 本につき，その炭素の OH を C=O に
    （2 本切れた炭素はギ酸 HCOOH）．戻り値: (消費 HIO4 = 切断数, 生成物 SMILES のリスト)"""
    m = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(smi)))
    def oh_of(c):
        for n in c.GetNeighbors():
            if n.GetSymbol() == 'O' and n.GetDegree() == 2 and sum(1 for x in n.GetNeighbors() if x.GetSymbol() == 'H') == 1:
                return n
        return None
    cuts = []
    for b in m.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetSymbol() == 'C' and a2.GetSymbol() == 'C' and oh_of(a1) and oh_of(a2):
            cuts.append((a1.GetIdx(), a2.GetIdx()))
    count = {}
    for i, j in cuts:
        m.RemoveBond(i, j)
        count[i] = count.get(i, 0) + 1; count[j] = count.get(j, 0) + 1
    # each cut carbon: lost C neighbor → add H? Instead: carbon with k cuts: OH->C=O (remove H on O, set double), and if k==2 add new OH
    for c, k in count.items():
        at = m.GetAtomWithIdx(c)
        o = oh_of(at)
        h = [n.GetIdx() for n in o.GetNeighbors() if n.GetSymbol() == 'H'][0]
        m.RemoveBond(o.GetIdx(), h)
        m.GetBondBetweenAtoms(c, o.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
        # 失った H（O 上）は捨て，C には失った C の代わりに H/OH をつける
        if k == 2:
            oid = m.AddAtom(Chem.Atom(8)); m.AddBond(c, oid, Chem.BondType.SINGLE)
        at.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    # remove orphan H atoms (from O-H removed)
    orphans = sorted([a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'H' and a.GetDegree() == 0], reverse=True)
    for i in orphans: m.RemoveAtom(i)
    mol = m.GetMol()
    for a in mol.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    frags = sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mol, asMols=True))
    return len(cuts), frags
import re
from rdkit.Chem import rdMolDescriptors as D
it = load('POL-SUG-022'); n = check_all_equations(it); check_subitems(it); assert n == 4
A = 'CO[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O'      # methyl alpha-D-glucopyranoside
B = 'COC1[C@H](O)[C@@H](O)[C@H](O1)[C@H](O)CO'        # methyl D-glucofuranoside
# 立体の確認：A を開環（環 O–C1 切断，OMe を除く）すると D-グルコース鎖状と同じ骨格
assert [c for _, c in Chem.FindMolChiralCenters(Chem.MolFromSmiles(A))] == ['S', 'R', 'S', 'S', 'R']
def flip(s, idx):
    t = list(re.finditer(r'\[C@@?H\]', s))[idx]
    return s[:t.start()] + ('[C@H]' if t.group(0) == '[C@@H]' else '[C@@H]') + s[t.end():]
MAN, GAL, BGLC = flip(A, 4), flip(A, 2), flip(A, 0)   # トークン順 C1,C5,C4,C3,C2
# マンノース(2位S)・ガラクトース(4位R) の CIP 確認
assert Chem.FindMolChiralCenters(Chem.MolFromSmiles(MAN))[-1][1] == 'S'
assert Chem.FindMolChiralCenters(Chem.MolFromSmiles(GAL))[2][1] == 'R'
def summary(s):
    k, fr = periodate(s)
    return k, fr, fr.count('O=CO'), fr.count('C=O')
kA, frA, fa, ha = summary(A); kB, frB, fb, hb = summary(B)
print('A', kA, frA); print('B', kB, frB)
assert (kA, fa, ha) == (2, 1, 0) and (kB, fb, hb) == (2, 0, 1)
assert sorted(D.CalcMolFormula(Chem.MolFromSmiles(f)) for f in frA) == ['C6H10O5', 'CH2O2']
assert sorted(D.CalcMolFormula(Chem.MolFromSmiles(f)) for f in frB) == ['C6H8O5', 'CH2O']
# (2)
M = 7*12 + 14*1.0 + 6*16; ntot = 1.94 / M; nA = 6.0e-3; nB = ntot - nA
print('(2)', ntot, nB * M, 2 * ntot, nB)
assert abs(ntot - 0.0100) < 1e-6 and sig(nB * M, 2) == 0.78 and sig(2*ntot, 2) == 0.020 and sig(nB, 2) == 0.0040
# (3) メチルキシロシド
XF = 'COC1[C@H](O)[C@@H](O)[C@H](O1)CO'; XP = 'COC1[C@H](O)[C@@H](O)[C@H](O)CO1'
Mx = 6*12 + 12*1.0 + 5*16; nx = 0.82 / Mx
assert D.CalcMolFormula(Chem.MolFromSmiles(XF)) == 'C6H12O5' and abs(nx - 5.0e-3) < 1e-6
kf, frf, ff, _ = summary(XF); kp, frp, fp, _ = summary(XP)
print('xyl-f', kf, frf, 'xyl-p', kp, frp)
assert (kf, ff) == (1, 0) and (kp, fp) == (2, 1)
# (4) 生成物の同一性
prod = {n: [f for f in periodate(s)[1] if f != 'O=CO'][0] for n, s in [('A', A), ('Man', MAN), ('Gal', GAL), ('bGlc', BGLC)]}
print(prod)
assert prod['A'] == prod['Man'] == prod['Gal'] != prod['bGlc']
for s_ in ['B $0.78\\,\\mathrm{g}$', '$\\ce{HIO4}$ $2.0\\times10^{-2}\\,\\mathrm{mol}$', '(3) 五員環']: in_text(s_, it)
print('POL-SUG-022 OK')
