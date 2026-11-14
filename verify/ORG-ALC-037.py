# 検証スクリプト（東工大 g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/ORG-ALC-037/item.json'
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
from collections import defaultdict
def alkanes(n):
    cur = {'C'}
    for _ in range(n-1):
        nxt = set()
        for s in cur:
            m = Chem.MolFromSmiles(s)
            for a in m.GetAtoms():
                if a.GetDegree() < 4:
                    rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(6)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE); nxt.add(Chem.MolToSmiles(rw))
        cur = nxt
    return cur
alc = set()
for s in alkanes(6):
    m = Chem.MolFromSmiles(s)
    for a in m.GetAtoms():
        if a.GetDegree() >= 4: continue
        rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(8)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE); alc.add(Chem.MolToSmiles(rw))
assert len(alc) == 17
def can(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s))
def cls(s):
    m = Chem.MolFromSmiles(s); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]; return o.GetNeighbors()[0].GetDegree() - 1
def iod(s): return Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CH3][CX4;H1][OX2H]'))
def nch(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
def alkenes(C):
    m = Chem.MolFromSmiles(C); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]; ca = o.GetNeighbors()[0]; out = set()
    for nb in ca.GetNeighbors():
        if nb.GetSymbol() == 'C' and nb.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m); rw.RemoveAtom(o.GetIdx()); i, j = ca.GetIdx(), nb.GetIdx(); i -= (i > o.GetIdx()); j -= (j > o.GetIdx())
            rw.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.DOUBLE); mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return out
def ct(s): return any(si.type == Chem.StereoType.Bond_Double for si in Chem.FindPotentialStereo(Chem.MolFromSmiles(s)))
def skel(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s))
    for b in m.GetBonds(): b.SetBondType(Chem.BondType.SINGLE)
    return can(Chem.MolToSmiles(m))
def key(s): return (cls(s) == 1, iod(s), nch(s) > 0, any(ct(x) for x in alkenes(s)))
g = defaultdict(set)
for s in alc: g[key(s)].add(s)
names = {can('CCCCC(C)O'): '2-hexanol', can('CCC(C)C(C)O'): '3-methyl-2-pentanol', can('CC(C)CC(C)O'): '4-methyl-2-pentanol'}
c2 = g[(False, True, True, True)]; print(c2); assert c2 == set(names)
# (3)
sk = {s: frozenset(skel(x) for x in alkenes(s)) for s in c2}; assert len(set(sk.values())) == 3 and all(len(v) == 1 for v in sk.values())
assert len({nch(s) for s in c2}) == 2      # (b) は不十分
def ketone(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s)); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    m.GetBondBetweenAtoms(o.GetIdx(), o.GetNeighbors()[0].GetIdx()).SetBondType(Chem.BondType.DOUBLE)
    mm = m.GetMol(); mm.GetAtomWithIdx(o.GetIdx()).SetNoImplicit(False); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
kc = [nch(ketone(s)) > 0 for s in c2]; print(kc); assert sorted(kc) == [False, False, True]   # (c) は不十分
assert len({len(alkenes(s)) for s in c2}) == 1      # (d) は不十分
# (4)
uniq = sorted(next(iter(v)) for v in g.values() if len(v) == 1); print(uniq)
assert set(uniq) == {can('CC(O)C(C)(C)C'), can('CCC(C)(O)CC')}
# (5)
g5 = g[(True, False, False, False)]; print(g5); assert len(g5) == 5
none = [s for s in g5 if len(alkenes(s)) == 0]; assert none == [can('CCC(C)(C)CO')]
assert sum(len(v) for v in g.values()) == 17

print('ALL OK')
