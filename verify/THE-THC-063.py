ITEM = '/home/claude/chemx/x1/out/items/THE-THC-063/item.json'
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
h = {'1-butene':127,'cis':120,'trans':116,'bd':237}
assert sorted(['1-butene','cis','trans'], key=lambda k:h[k]) == ['trans','cis','1-butene']
chk('isom', -127-(-116), -11, 2)
chk('est', 2*127, 254, 3); chk('stab', 254-237, 17, 2)
assert 239.7 > 231.5
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
for s,f in [('C=CCC','C4H8'),('C/C=C\\C','C4H8'),('C=CC=C','C4H6'),('C1=CCCC=C1','C6H8'),('C1=CCC=CC1','C6H8')]:
    assert F(Chem.MolFromSmiles(s))==f, s
for e in ['C4H8 + H2 -> C4H10','C4H6 + 2H2 -> C4H10','C6H8 + 2H2 -> C6H12']:
    assert check_equation(e) is None, e
print('formulas OK')

print('ALL OK')
