# 検証スクリプト（東工大 g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/ORG-STR-061/item.json'
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
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def can(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s))
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
def nCH3(s):
    m = Chem.MolFromSmiles(s)
    return sum(1 for a in m.GetAtoms() if a.GetSymbol()=='C' and a.GetTotalNumHs()==3)
def alkyls(n):
    out = set()
    for s in alkanes(n):
        m = Chem.MolFromSmiles(s)
        for a in m.GetAtoms():
            rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(0)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
            out.add(Chem.MolToSmiles(rw))
    return out
def lab(R, n): return R.replace('[*]','*').replace('*', '[*:%d]' % n)
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
h6 = alkanes(6); assert len(h6) == 5
c4 = [s for s in h6 if nCH3(s) == 4]; print(c4); assert len(c4) == 2
# C5H10
alk5 = set()
for s in alkanes(5):
    m = Chem.MolFromSmiles(s)
    for b in m.GetBonds():
        rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
        try: Chem.SanitizeMol(rw); alk5.add(Chem.MolToSmiles(rw))
        except Exception: pass
cyc5 = {can(x) for x in ['C1CCCC1','CC1CCC1','CC1(C)CC1','CC1CC1C','CCC1CC1']}
all5 = alk5 | cyc5; assert len(alk5) == 5 and len(all5) == 10
assert all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C5H10' for s in all5)
zero = [s for s in all5 if nCH3(s) == 0]; print(zero); assert zero == [can('C1CCCC1')]
# C4H8O2 エステル  R-C(=O)-O-R'
est = set()
for k in range(0, 3):
    Rs = ['[H]*'] if k == 0 else alkyls(k)
    for R in Rs:
        for R2 in alkyls(3-k):
            m = Chem.molzip(Chem.MolFromSmiles('[*:1]C(=O)O[*:2].' + lab(R,1) + '.' + lab(R2,2)))
            est.add(Chem.MolToSmiles(Chem.RemoveHs(m)))
assert len(est) == 4 and all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H8O2' for s in est)
one = [s for s in est if nCH3(s) == 1]; print(one); assert one == [can('O=COCCC')]
# C5H12O エーテル
eth = set()
for k in (1, 2):
    for R in alkyls(k):
        for R2 in alkyls(5-k):
            m = Chem.molzip(Chem.MolFromSmiles('[*:1]O[*:2].' + lab(R,1) + '.' + lab(R2,2)))
            eth.add(Chem.MolToSmiles(m))
assert len(eth) == 6
mx = max(nCH3(s) for s in eth); top = [s for s in eth if nCH3(s) == mx]; print(top, mx)
assert mx == 4 and top == [can('COC(C)(C)C')]
h8 = alkanes(8); assert len(h8) == 18
mx8 = max(nCH3(s) for s in h8); top8 = [s for s in h8 if nCH3(s) == mx8]; print(top8, mx8)
assert mx8 == 6 and top8 == [can('CC(C)(C)C(C)(C)C')]

print('ALL OK')
