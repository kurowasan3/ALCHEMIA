ITEM = '/home/claude/chemx/x2/outh/items/ORG-CBX-074/item.json'
# 検証スクリプト（x2 g03）  python3 -I この.py
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
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS
# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import combinations
# (1) 候補の網羅：C5 骨格の末端に COOH 2 個，炭素 2 個と結合した炭素に C=O
res = set()
for s in ['CCCCC', 'CC(C)CC', 'CC(C)(C)C']:
    m = Chem.MolFromSmiles(s); n = m.GetNumAtoms()
    for cooh in combinations(range(n), 2):
        if any(m.GetAtomWithIdx(i).GetDegree() != 1 for i in cooh): continue
        for k in range(n):
            if k in cooh or m.GetAtomWithIdx(k).GetDegree() != 2: continue
            rw = Chem.RWMol(m)
            for i in cooh:
                o1 = rw.AddAtom(Chem.Atom(8)); o2 = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(i, o1, Chem.BondType.DOUBLE); rw.AddBond(i, o2, Chem.BondType.SINGLE)
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(k, o, Chem.BondType.DOUBLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); res.add(Chem.MolToSmiles(mm))
cand = {'a': 'CC(=O)C(C(=O)O)C(=O)O', 'b': 'OC(=O)C(C)C(=O)C(=O)O', 'c': 'OC(=O)CC(=O)CC(=O)O', 'd': 'OC(=O)C(=O)CCC(=O)O'}
can = {k: Chem.CanonSmiles(v) for k, v in cand.items()}
assert set(can.values()) == res and len(res) == 4, res
for k, v in can.items(): assert CalcMolFormula(Chem.MolFromSmiles(v)) == 'C5H6O5'
print('候補4つ OK')
# (2) モデルによる脱炭酸を機械的に適用
beta = Chem.MolFromSmarts('[CX3](=O)([#6])[CX4][CX3](=O)[OX2H1]')
malon = Chem.MolFromSmarts('[CX4]([CX3](=O)[OX2H1])[CX3](=O)[OX2H1]')
def decarb(smi):
    m = Chem.MolFromSmiles(smi); n = 0
    while True:
        hit = None
        for patt, ci, oi in [(beta, 4, 6), (malon, 1, 3)]:
            mt = m.GetSubstructMatches(patt)
            if mt: hit = (mt[0][ci], mt[0][ci+1], mt[0][ci+2]); break
        if not hit: return n, Chem.MolToSmiles(m)
        rw = Chem.RWMol(m)
        for idx in sorted(hit, reverse=True): rw.RemoveAtom(idx)
        m = rw.GetMol(); Chem.SanitizeMol(m); m = Chem.MolFromSmiles(Chem.MolToSmiles(m)); n += 1
exp = {'a': (2, 'CC(C)=O'), 'b': (1, 'CCC(=O)C(=O)O'), 'c': (2, 'CC(C)=O'), 'd': (0, can['d'])}
for k, v in can.items():
    r = decarb(v); assert r == (exp[k][0], Chem.CanonSmiles(exp[k][1])), (k, r); print(k, r)
# (3)(4) 不斉炭素
def nst(smi): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(smi), includeUnassigned=True, useLegacyImplementation=False))
assert [nst(can[k]) for k in 'abcd'] == [0, 1, 0, 0]
red = {'a': 'CC(O)C(C(=O)O)C(=O)O', 'c': 'OC(=O)CC(O)CC(=O)O', 'b': 'OC(=O)C(C)C(O)C(=O)O'}
assert nst(red['a']) == 1 and nst(red['c']) == 0 and nst(red['b']) == 2
assert nst('CCC(=O)C(=O)O') == 0
# b の還元物の立体異性体を列挙（メソなし→4）
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
iso = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(red['b']))}
assert len(iso) == 4, iso
print('立体異性体', len(iso))
print('ALL OK')
