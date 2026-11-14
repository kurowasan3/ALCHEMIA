ITEM = '/home/claude/chemx/x1/outh/items/ORG-HYD-053/item.json'
# 検証スクリプト（x1 g05）  python3 -I この.py
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS, it['field']
for f in it['format'].split(','):
    assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
import os
d = os.path.dirname(ITEM)
for k in ('problem_fig','solution_fig'):
    if it[k]: assert os.path.exists(os.path.join(d, it[k])), it[k]

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
def place(smi, order):
    m = Chem.MolFromSmiles(smi); res = set()
    for b in m.GetBonds():
        rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(order)
        try: Chem.SanitizeMol(rw)
        except Exception: continue
        res.add(Chem.MolToSmiles(rw))
    return res
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
alk = ['CCCCCC','CC(C)CCC','CCC(C)CC','CC(C)C(C)C','CC(C)(C)CC']
cyc = ['C1CCCCC1','CC1CCCC1']
ALL = sorted(set().union(*[place(s, Chem.BondType.TRIPLE) for s in alk]) | set().union(*[place(s, Chem.BondType.DOUBLE) for s in cyc]))
assert len(ALL) == 12 and all(F(Chem.MolFromSmiles(s)) == 'C6H10' for s in ALL)
print('候補 12 種類 OK')
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
def multi(m):
    return [b for b in m.GetBonds() if b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)][0]
def stereo(s):
    m = Chem.MolFromSmiles(s)
    return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
def cistrans_possible(s):
    m = Chem.MolFromSmiles(s)
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
    return len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)}) > 1 and stereo(s) == 0
def hyd(s, full):
    m = Chem.RWMol(Chem.MolFromSmiles(s)); b = multi(m)
    b.SetBondType(Chem.BondType.SINGLE if full or b.GetBondType()==Chem.BondType.DOUBLE else Chem.BondType.DOUBLE)
    mm = m.GetMol(); mm.UpdatePropertyCache(strict=False); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
def branched(s):
    m = Chem.MolFromSmiles(s); return any(sum(1 for n in a.GetNeighbors() if n.GetSymbol()=='C') > 2 for a in m.GetAtoms() if a.GetSymbol()=='C') or m.GetRingInfo().NumRings() > 0
def terminal_alkyne(s): return Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CH]#C'))
from rdkit.Chem import AllChem
_rx = AllChem.ReactionFromSmarts('[CX3;H1,H2:1]=[O:2]>>[C:1](=[O:2])O')
def cleave(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s)); b = multi(m)
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); m.RemoveBond(i, j)
    for a in (i, j):
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(a, o, Chem.BondType.DOUBLE)
    m = m.GetMol(); m.UpdatePropertyCache(strict=False); Chem.SanitizeMol(m)
    out = []
    for sm in Chem.MolToSmiles(m).split('.'):
        fm = Chem.MolFromSmiles(sm)
        while fm.HasSubstructMatch(Chem.MolFromSmarts('[CX3;H1,H2]=O')):
            fm = _rx.RunReactants((fm,))[0][0]; Chem.SanitizeMol(fm)
        c = Chem.MolToSmiles(fm)
        if c in (can('OC(=O)O'), can('O=CO')): c = 'O=C=O'
        out.append(c)
    return sorted(out)
def cleave_tri(s):
    # alkyne: each sp carbon -> COOH (terminal CH -> CO2)
    m = Chem.RWMol(Chem.MolFromSmiles(s)); b = multi(m); i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); m.RemoveBond(i, j)
    for a in (i, j):
        o1 = m.AddAtom(Chem.Atom(8)); m.AddBond(a, o1, Chem.BondType.DOUBLE)
        o2 = m.AddAtom(Chem.Atom(8)); m.AddBond(a, o2, Chem.BondType.SINGLE)
    m = m.GetMol(); m.UpdatePropertyCache(strict=False); Chem.SanitizeMol(m)
    out = []
    for sm in Chem.MolToSmiles(m).split('.'):
        out.append('O=C=O' if can(sm) in (can('O=C(O)O'), can('O=CO')) else can(sm))
    return sorted(out)
def ox(s):
    return cleave_tri(s) if multi(Chem.MolFromSmiles(s)).GetBondType()==Chem.BondType.TRIPLE else cleave(s)
def methylketone(s): return Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)[#6]'))
def hydration(s):
    m0 = Chem.MolFromSmiles(s); b = multi(m0)
    if b.GetBondType() != Chem.BondType.TRIPLE: return None
    out = set()
    for a in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()):
        m = Chem.RWMol(m0); m.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx()).SetBondType(Chem.BondType.SINGLE)
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(a, o, Chem.BondType.DOUBLE); mm = m.GetMol(); mm.UpdatePropertyCache(strict=False); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return out
cyclic = lambda s: Chem.MolFromSmiles(s).GetRingInfo().NumRings() > 0
P = [s for s in ALL if terminal_alkyne(s) and stereo(s) > 0]
Q = [s for s in ALL if len(set(p for p in ox(s) if p!='O=C=O'))==1 and 'O=C=O' not in ox(s) and not branched(ox(s)[0]) and methylketone(ox(s)[0])]
R = [s for s in ALL if not terminal_alkyne(s) and cistrans_possible(hyd(s, False)) and not branched(hyd(s, True)) and hydration(s) is not None and len(hydration(s))==1]
S = [s for s in ALL if not terminal_alkyne(s) and stereo(s) > 0 and stereo(hyd(s, True)) == 0]
print('P', P, 'Q', Q, 'R', R, 'S', S)
assert P == [can('C#CC(C)CC')] and Q == [can('CC1=CCCC1')] and R == [can('CCC#CCC')], (P, Q, R)
assert S == [can('CC1C=CCC1')], S
qo = ox(Q[0]); assert qo == [can('CC(=O)CCCC(=O)O')], qo
so = ox(S[0]); print('S ox', so); assert so == [can('CC(CCC(=O)O)C(=O)O')] and stereo(so[0]) == 1
# 2-hexyne gives two ketones
assert len(hydration(can('CC#CCCC'))) == 2
# without the "cyclic" constraint, check which compounds satisfy each property, to report
print('stereocenter compounds:', [s for s in ALL if stereo(s) > 0])
print('一意性 OK')

print('ALL OK')
