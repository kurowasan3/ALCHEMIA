ITEM = '/home/claude/chemx/x1/out/items/ORG-CBX-066/item.json'
# 検証スクリプト（x1 g02）  python3 -I この.py
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
FIELDS = ['原子の構造と電子配置','元素の周期表と周期性','化学結合と結晶','物質量と化学反応式','酸と塩基','酸化還元反応','電池と電気分解','気体の性質','溶液の性質','化学反応とエネルギー','反応速度と化学平衡','電離平衡と溶解平衡','ハロゲン・貴ガス・水素','酸素・硫黄','窒素・リン','炭素・ケイ素','アルカリ金属・アルカリ土類金属','アルミニウム・亜鉛・スズ・鉛','遷移元素（鉄・銅・銀・クロム・マンガン）','金属イオンの沈殿・錯イオン・系統分析','気体の製法と性質','有機化合物の構造と元素分析','脂肪族炭化水素','アルコール・アルデヒド・ケトン','カルボン酸・エステル・油脂','芳香族化合物','有機化合物の検出反応と分離','高分子の基礎と合成樹脂','合成繊維','ゴム・機能性高分子','糖類','アミノ酸・タンパク質・核酸']
assert it['field'] in FIELDS, it['field']

# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import FindMolChiralCenters
assert check_equation('C6H5CH=CHCOOCH(C2H5)2 + H2O -> C6H5CH=CHCOOH + (C2H5)2CHOH') is None
X = Chem.MolFromSmiles('O=C(/C=C/c1ccccc1)OC(CC)CC')
assert CalcMolFormula(X) == 'C14H18O2'
M = 14*12 + 18*1.0 + 2*16; assert M == 218
n = 43.6/M
chk('CO2', n*14*44, 123.2, 4); chk('H2O', n*9*18, 32.4, 3)
# C9H8O2 の候補とシス-トランス
cands = ['OC(=O)C=Cc1ccccc1','OC(=O)C(=C)c1ccccc1','OC(=O)c1ccccc1C=C','OC(=O)c1cccc(C=C)c1','OC(=O)c1ccc(C=C)cc1']
from rdkit.Chem import EnumerateStereoisomers as ES
ct = []
for s in cands:
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C9H8O2'
    ct.append(len(list(ES.EnumerateStereoisomers(m))))
assert ct == [2,1,1,1,1], ct
# 第二級アルコール C5H12O
sec = {'2-ペンタノール':'CC(O)CCC','3-ペンタノール':'CCC(O)CC','3-メチル-2-ブタノール':'CC(O)C(C)C'}
patt = Chem.MolFromSmarts('[CH3][CH1]([OH])')
for k, s in sec.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C5H12O'
    has = m.HasSubstructMatch(patt)
    nch = len(FindMolChiralCenters(m, includeUnassigned=True))
    print(k, 'CH3CH(OH)', has, '不斉', nch)
    assert has == (k != '3-ペンタノール') and nch == (0 if k == '3-ペンタノール' else 1)
for s in ['O=C(C=Cc1ccccc1)OC(CC)CC', 'O=C(C=Cc1ccccc1)OC(C)CCC', 'O=C(C=Cc1ccccc1)OC(C)C(C)C']:
    print(s, len(list(ES.EnumerateStereoisomers(Chem.MolFromSmiles(s)))))
assert [len(list(ES.EnumerateStereoisomers(Chem.MolFromSmiles(s)))) for s in ['O=C(C=Cc1ccccc1)OC(CC)CC', 'O=C(C=Cc1ccccc1)OC(C)CCC', 'O=C(C=Cc1ccccc1)OC(C)C(C)C']] == [2,4,4]

print('ALL OK')
