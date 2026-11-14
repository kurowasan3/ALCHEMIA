# 検証スクリプト（東工大 g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/ORG-CBX-064/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
FIELDS = ['原子の構造と電子配置','元素の周期表と周期性','化学結合と結晶','物質量と化学反応式','酸と塩基','酸化還元反応','電池と電気分解','気体の性質','溶液の性質','化学反応とエネルギー','反応速度と化学平衡','電離平衡と溶解平衡','ハロゲン・貴ガス・水素','酸素・硫黄','窒素・リン','炭素・ケイ素','アルカリ金属・アルカリ土類金属','アルミニウム・亜鉛・スズ・鉛','遷移元素（鉄・銅・銀・クロム・マンガン）','金属イオンの沈殿・錯イオン・系統分析','気体の製法と性質','有機化合物の構造と元素分析','脂肪族炭化水素','アルコール・アルデヒド・ケトン','カルボン酸・エステル・油脂','芳香族化合物','有機化合物の検出反応と分離','高分子の基礎と合成樹脂','合成繊維','ゴム・機能性高分子','糖類','アミノ酸・タンパク質・核酸']
assert it['field'] in FIELDS, it['field']

# ---- 個別の検証 ----
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def alkanes(n):
    cur = {'C'}
    for _ in range(n-1):
        nxt = set()
        for s in cur:
            m = Chem.MolFromSmiles(s)
            for a in m.GetAtoms():
                if a.GetDegree() < 4:
                    rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(6)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
                    nxt.add(Chem.MolToSmiles(rw))
        cur = nxt
    return cur
def subst(n, grp, k):
    cur = alkanes(n)
    for _ in range(k):
        nxt = set()
        for s in cur:
            m = Chem.MolFromSmiles(s); N = m.GetNumAtoms()
            for a in m.GetAtoms():
                if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0:
                    rw = Chem.RWMol(Chem.CombineMols(m, Chem.MolFromSmiles(grp)))
                    rw.AddBond(a.GetIdx(), N, Chem.BondType.SINGLE)
                    mm = rw.GetMol(); Chem.SanitizeMol(mm); nxt.add(Chem.MolToSmiles(mm))
        cur = nxt
    return cur
def nch(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
acids = {d: {s for s in subst(d-2, 'C(=O)O', 2) if nch(s) > 0} for d in (4, 5, 6, 7)}
alcs = {n: subst(n, 'O', 1) for n in range(1, 7)}
print({d: len(v) for d, v in acids.items()})
assert len(acids[4]) == 0 and len(acids[5]) == 1 and len(acids[6]) == 3
assert all(nch(s) == 0 for s in alcs[3]) and sum(nch(s) > 0 for s in alcs[4]) == 1
# A の組み立て：D の 2 つの O-H の O に B, C をつなぐ
def esters(D, B, C):
    out = set(); m = Chem.MolFromSmiles(D)
    ohs = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'O' and a.GetTotalNumHs() == 1]
    assert len(ohs) == 2
    for o1, o2 in ((ohs[0], ohs[1]), (ohs[1], ohs[0])):
        rw = Chem.RWMol(m)
        for o, R in ((o1, B), (o2, C)):
            r = Chem.MolFromSmiles(R); N = rw.GetNumAtoms()
            ro = [a.GetIdx() for a in r.GetAtoms() if a.GetSymbol() == 'O'][0]
            # アルコールの O を除いて C をつなぐ
            rr = Chem.RWMol(r); cidx = [n.GetIdx() for n in r.GetAtomWithIdx(ro).GetNeighbors()][0]
            rr.RemoveAtom(ro); cidx = cidx if cidx < ro else cidx - 1
            N = rw.GetNumAtoms(); rw = Chem.RWMol(Chem.CombineMols(rw, rr.GetMol())); rw.AddBond(o, N + cidx, Chem.BondType.SINGLE)
        mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return out
combos = []
for d in range(4, 12):
    for c in range(1, 12):
        b = 11 - c - d
        if b < 1: continue
        Ds = acids.get(d, set()) if d <= 7 else set()
        Cs = [s for s in alcs.get(c, []) if nch(s) > 0]
        Bs = [s for s in alcs.get(b, []) if nch(s) == 0]
        if Ds and Cs and Bs: combos.append((b, c, d))
print(combos); assert sorted(combos) == sorted([(2,4,5),(1,5,5),(1,4,6)])
allA = {}
for b, c, d in combos:
    for D in acids[d]:
        for C in [s for s in alcs[c] if nch(s) > 0]:
            for B in [s for s in alcs[b] if nch(s) == 0]:
                for A in esters(D, B, C): allA[A] = (B, C, D)
assert all(CalcMolFormula(Chem.MolFromSmiles(a)) == 'C11H20O4' for a in allA)
print(len(allA)); assert len(allA) == 13
def succinic(D):  # 2 つの COOH が隣り合う炭素に結合
    return Chem.MolFromSmiles(D).HasSubstructMatch(Chem.MolFromSmarts('OC(=O)[CX4][CX4]C(=O)O'))
def secondary(C):
    return Chem.MolFromSmiles(C).HasSubstructMatch(Chem.MolFromSmarts('[CX4;H1]([#6])([#6])[OX2H]'))
def iodoform(B):
    return Chem.MolFromSmiles(B).HasSubstructMatch(Chem.MolFromSmarts('[CH3][CX4;H1,H2][OX2H]'))
A3 = {a: v for a, v in allA.items() if succinic(v[2]) and secondary(v[1]) and not iodoform(v[0])}
print(len(A3)); assert len(A3) == 7
def alkenes(C):
    m = Chem.MolFromSmiles(C); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    ca = o.GetNeighbors()[0]; out = set()
    for nb in ca.GetNeighbors():
        if nb.GetSymbol() == 'C' and nb.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); rw.RemoveAtom(o.GetIdx())
            i, j = ca.GetIdx(), nb.GetIdx(); i -= (i > o.GetIdx()); j -= (j > o.GetIdx())
            rw.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.DOUBLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return out
def has_ct(s):
    m = Chem.MolFromSmiles(s)
    from rdkit.Chem import FindPotentialStereo
    return any(si.type == Chem.StereoType.Bond_Double for si in Chem.FindPotentialStereo(m))
A4 = {a: v for a, v in A3.items() if not any(has_ct(x) for x in alkenes(v[1]))}
print(A4); assert len(A4) == 2
ii = Chem.MolToSmiles(Chem.MolFromSmiles('CC(C)C(C)OC(=O)CC(C)C(=O)OC'))
i_ = Chem.MolToSmiles(Chem.MolFromSmiles('CC(C)C(C)OC(=O)C(C)CC(=O)OC'))
assert set(A4) == {i_, ii}
# (5) メチルエステルを外した E：COOH の結合炭素が不斉か
Eii = 'CC(C)C(C)OC(=O)CC(C)C(=O)O'; Ei = 'CC(C)C(C)OC(=O)C(C)CC(=O)O'
def cooh_on_stereo(s):
    m = Chem.MolFromSmiles(s); cc = dict(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
    patt = Chem.MolFromSmarts('[CX4]C(=O)[OX2H]')
    return any(t[0] in cc for t in m.GetSubstructMatches(patt))
assert cooh_on_stereo(Eii) and not cooh_on_stereo(Ei)

print('ALL OK')
