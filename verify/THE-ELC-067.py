# 検証スクリプト（東工大 g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/THE-ELC-067/item.json'
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
F = 96500
ne = 0.500*1158/F; assert abs(ne-6.00e-3) < 1e-9
nAg = 2.00e-3; nCu = (ne-nAg)/2; assert nCu < 3.00e-3   # Cu2+ は残る → H2 は出ない
m = 108*nAg + 63.5*nCu; print(m); assert sig(m,3) == 0.343
V = ne/4*22.4e3; assert sig(V,3) == 33.6
t = nAg*F/0.500; assert sig(t,3) == 386
r = (63.5/2)/108; print(r); assert sig(r,2) == 0.29
for e in ['Ag+ + e- -> Ag','Cu^2+ + 2e- -> Cu','2H2O -> O2 + 4H+ + 4e-']:
    assert check_equation(e) is None

print('ALL OK')
