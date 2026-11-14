# 検証スクリプト（x1 g07）  python3 -I THE-ION-109.py
ITEM = '/home/claude/chemx/x1/outh/items/THE-ION-109/item.json'
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


import sympy as sp
e=65000/(2.3*8.31)*(1/298-1/333); assert round(e,2)==1.20
K2=1.8e-10*10**e; assert sig(K2,2)==2.8e-9
s60=sqrt(K2); s25=sqrt(1.8e-10); assert sig((s60-s25)*143.5*1000,2)==5.7
assert 1e-4*2e-5 < K2 and 1e-4*2e-5 > 1.8e-10
x=sp.symbols('x'); sols=[float(r) for r in sp.solve((1e-4-x)*(2e-5-x)-1.8e-10,x) if 0<float(r)<2e-5]
assert len(sols)==1 and round(sols[0]/2e-5*100)==89
assert sig(10**(3900/(2.3*8.31)*(1/298-1/333)),2)==1.2
assert round(10**0.20,2)==1.58 and round(sqrt(1.8),2)==1.34 and round(sqrt(28.5),2)==5.34 and round(sqrt(71.2),2)==8.44 and round(10**0.072,2)==1.18
# 格子エンタルピーのサイクル
U=284+735+243/2-349+127; assert abs(U-918.5)<0.1 and abs(U-513-340-65.5)<0.1

print('ALL OK')
