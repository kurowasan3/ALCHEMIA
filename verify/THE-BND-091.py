ITEM = '/home/claude/chemx/x1/outh/items/THE-BND-091/item.json'
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

from itertools import product
from collections import Counter
def cluster(K):
    return [p for p in product(range(-2*K,2*K+1),repeat=3) if sum(p)%2==0 and shell(p) <= K]
def shell(p):
    # 立方八面体：正方形の面 |x_i|=k，正三角形の面 |x|+|y|+|z|=2k（最近接距離 √2 の座標）
    return max(max(abs(t) for t in p), (sum(abs(t) for t in p)+1)//2)
for K, N, S in [(1,13,12),(2,55,42),(3,147,92)]:
    c = cluster(K); assert len(c) == N, (K, len(c))
    assert sum(1 for p in c if shell(p) == K) == S == 10*K*K+2
def cns(K):
    c = set(cluster(K)); out = Counter()
    for p in c:
        if shell(p) != K: continue
        n = sum(1 for q in c if sum((a-b)**2 for a,b in zip(p,q)) == 2)
        out[n] += 1
    return out
c2 = cns(2); print('55:', c2); assert c2 == Counter({7:24, 5:12, 8:6})
c3 = cns(3); print('147:', c3); assert c3[5] + c3[7] == 60 and c3[5] == 12 and c3[7] == 48
chk('N', 5.50e-3*6.02e23/6.02e19, 55.0, 3)
chk('surf', 42/55*100, 76, 2)
chk('ratio', 36*147/(55*60), 1.6, 2)

print('ALL OK')
