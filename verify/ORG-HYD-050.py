ITEM = '/home/claude/chemx/x1/out/items/ORG-HYD-050/item.json'
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
from rdkit.Chem import EnumerateStereoisomers as ES
from itertools import combinations_with_replacement
halves = ['C', 'CC', 'CCC']      # =CH2, =CHCH3, =CHC2H5 の炭素鎖
seen = {}
for a, b in combinations_with_replacement(halves, 2):
    m = Chem.MolFromSmiles(('C'*(len(a)-1)) + 'C=C' + ('C'*(len(b)-1)))
    seen[Chem.MolToSmiles(m)] = m
react = {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in ['C=CC', 'C=CCC']}
new = [k for k in seen if k not in react]
assert len(seen) == 6 and len(new) == 4, (seen.keys(), new)
tot = sum(len(list(ES.EnumerateStereoisomers(seen[k]))) for k in new)
assert tot == 7, tot
print('新しいアルケン', new, '立体を区別して', tot)
assert check_equation('2CH2=CHCH2CH3 <=> CH2=CH2 + CH3CH2CH=CHCH2CH3') is None
m = Chem.MolFromSmiles('C1=CCCCC1'); assert CalcMolFormula(m) == 'C6H10'
assert m.GetRingInfo().AtomRings()[0].__len__() == 6

print('ALL OK')
