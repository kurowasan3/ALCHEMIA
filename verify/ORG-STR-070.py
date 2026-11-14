# ORG-STR-070 検証：C9H12O 芳香族異性体の網羅，不斉炭素をもつもの，各条件での絞り込み
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D

def grow(smis, n):
    for _ in range(n):
        new = set()
        for s in smis:
            m = Chem.MolFromSmiles(s)
            for a in m.GetAtoms():
                if a.GetTotalNumHs() == 0: continue
                for z in (6, 8):
                    if z == 8 and a.GetAtomicNum() == 8: continue   # O-O は除く
                    rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(z))
                    rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
                    try:
                        mm = rw.GetMol(); Chem.SanitizeMol(mm); new.add(Chem.MolToSmiles(mm))
                    except Exception: pass
        smis = new
    return smis
# 不飽和度 4 = ベンゼン環のみ。残り C3 と O1 は飽和の鎖
allm = [s for s in grow({'c1ccccc1'}, 4) if D.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C9H12O']
print('C9H12O aromatic isomers:', len(allm))
def nchiral(s):
    return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
chir = sorted(s for s in allm if nchiral(s) > 0)
print('chiral:', len(chir), chir); assert len(chir) == 7

def has(s, smarts): return Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts(smarts))
info = {}
for s in chir:
    phenol = has(s, 'c[OX2H]'); alcohol = has(s, '[CX4][OX2H]'); ether = has(s, '[#6][OX2][#6]')
    prim = has(s, '[CH2][OX2H]'); sec = has(s, '[CH1]([#6])([#6])[OX2H]')
    iodo = has(s, '[CH3][CH1][OX2H]')  # CH3CH(OH)-
    # 環の一塩素置換体の数
    m = Chem.MolFromSmiles(s); ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    ringH = {ranks[a.GetIdx()] for a in m.GetAtoms() if a.GetIsAromatic() and a.GetTotalNumHs() > 0}
    # 酸化生成物
    rxn_p = Chem.AllChem.ReactionFromSmarts('[CH2:1][OH:2]>>[C:1](=O)[OH:2]') if False else None
    info[s] = dict(phenol=phenol, alcohol=alcohol, ether=ether, prim=prim, sec=sec, iodo=iodo, ringCl=len(ringH))
    print(s, info[s])
from rdkit.Chem import AllChem
def oxidize(s):
    m = Chem.MolFromSmiles(s)
    for sm in ('[CH2:1][OH1:2]>>[C:1](=O)[OH1:2]', '[CH1:1]([#6:3])([#6:4])[OH1:2]>>[C:1]([#6:3])([#6:4])=[O:2]'):
        r = AllChem.ReactionFromSmarts(sm); ps = r.RunReactants((m,))
        if ps:
            p = ps[0][0]; Chem.SanitizeMol(p); return Chem.MolToSmiles(p)
    return None
achiral_ox = 0
for s in chir:
    p = oxidize(s)
    print('oxidation', s, '->', p, (nchiral(p) if p else None))
    if p and nchiral(p) == 0: achiral_ox += 1
print('oxidizable with achiral product:', achiral_ox); assert achiral_ox == 5
# 条件による同定
P = [s for s in chir if info[s]['prim']]
Q = [s for s in chir if not info[s]['alcohol'] and not info[s]['phenol']]
R = [s for s in chir if info[s]['sec'] and info[s]['iodo'] and info[s]['ringCl'] == 3]
S = [s for s in chir if info[s]['sec'] and not info[s]['iodo']]
T = [s for s in chir if info[s]['sec'] and info[s]['iodo'] and info[s]['ringCl'] == 2]
print(P, Q, R, S, T)
assert P == ['CC(CO)c1ccccc1'] and Q == ['COC(C)c1ccccc1'] and R == ['CC(O)Cc1ccccc1']
assert S == ['CCC(O)c1ccccc1'] and T == ['Cc1ccc(C(C)O)cc1']
assert nchiral(oxidize(P[0])) == 1
# 脱水生成物（シス-トランスを区別）
def dehyd(s):
    m = Chem.MolFromSmiles(s); out = set()
    r = AllChem.ReactionFromSmarts('[C:1]([OH1:3])-[C!H0:2]>>[C:1]=[C:2].[O:3]')
    for ps in r.RunReactants((m,)):
        p = ps[0]
        try:
            Chem.SanitizeMol(p)
            from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
            for iso in EnumerateStereoisomers(p): out.add(Chem.MolToSmiles(iso))
        except Exception: pass
    return out
for name, s in (('R', R[0]), ('S', S[0])):
    pr = dehyd(s); print(name, len(pr), pr)
assert len(dehyd(R[0])) == 3 and len(dehyd(S[0])) == 2
for s in (P[0], Q[0], R[0], S[0], T[0]):
    assert D.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C9H12O'
# 燃焼の式
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['2C9H12O + 23O2 -> 18CO2 + 12H2O']:
    print(e, check_equation(e)); assert check_equation(e) is None
print('OK')
