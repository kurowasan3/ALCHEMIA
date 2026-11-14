# 検証スクリプト（x1 g07）  python3 -I THE-ELC-081.py
ITEM = '/home/claude/chemx/x1/out/items/THE-ELC-081/item.json'
import sys, json, re
from math import log10, floor, log, exp, sqrt
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
assert it['field'] in FIELDS, it['field']
assert set(it['format'].split(',')) <= {'qa','equation','compare','condition','pathway','calc','structure'}

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

neq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            r = check_equation(b)
            assert r is None, (b, r)
            neq += 1
print('反応式チェック:', neq, '本 OK')
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans, (nprob, nans)
print('小問', nprob, '答え', nans)

def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def near(a, b, rel=0.01):
    return abs(a - b) <= rel * abs(b)


F = 9.65e4
Q = 0.50*7720; ne = Q/F
assert near(ne, 0.0400, 1e-6)
m = ne/4*32
assert sig(m,2) == 0.32
nZn = 0.65/65; Q3 = 2*nZn*F
assert sig(Q3,2) == 1900
V = nZn*2/4*22.4*1000
assert sig(V,2) == 110
# 半反応の和で OH-, H2O が消えること
assert check_equation('2Zn + 4OH- + O2 + 2H2O -> 2ZnO + 2H2O + 4OH-') is None
assert check_equation('2Zn + O2 -> 2ZnO') is None

print('ALL OK')
