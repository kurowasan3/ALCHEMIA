# 検証スクリプト（グループ g06）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

FIELDS = set('原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split())

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    for f in it['format'].split(','):
        assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    assert it['answer'].count('$') % 2 == 0 and it['solution'].count('$') % 2 == 0 and it['problem'].count('$') % 2 == 0
    return it

def sig(x, n=2):
    from math import log10, floor
    if x == 0: return 0
    d = n - 1 - floor(log10(abs(x)))
    return round(x, d)

it = check_item('/home/claude/chem6/outh/items/THE-ION-076/item.json')
from sympy import symbols, solve, sqrt, Rational as Q, nsolve
h, c = symbols('h c', positive=True)
K1S = Q(1,10**8); K12 = Q(1,10**21); S = Q(1,10)
hc = lambda cv: (cv + sqrt(cv**2 + 4*K1S))/2
vals = [float(hc(Q(0))), float(hc(Q(1,10**4))), float(hc(Q(1,100)))]
print(vals); assert [sig(v) for v in vals] == [1.0e-4, 1.6e-4, 1.0e-2]
print('和の近似のずれ', 2e-4/vals[1]); assert round(2e-4/vals[1],2) == 1.24
# 厳密（HS- の電離・水の電離も含む）でもほぼ同じ
Kw = Q(1,10**14); K2 = Q(1,10**14)
he = nsolve(h - Q(1,10**4) - K1S/h - 2*K1S*K2/h**2 - Kw/h, h, 1.6e-4); print('厳密', he); assert abs(float(he)/vals[1]-1) < 1e-4
# (3)
c3 = solve(sp := (K1S/h**2 - Q(1,100)), h)[0]; c3v = c3 - K1S/c3; print('(3)', float(c3v)); assert sig(float(c3v)) == 9.9e-4
# (4)
ip = lambda hv: Q(1,100)*K12*S/hv**2
print(float(ip(hc(Q(1,10**4)))), float(ip(Q(2,10**4))))
Ksp = Q(3,10**17)
assert ip(hc(Q(1,10**4))) > Ksp and ip(Q(2,10**4)) < Ksp and sig(float(ip(hc(Q(1,10**4))))) == 3.8e-17
# (5)
hth = sqrt(Q(1,100)*K12*S/Ksp); cth = hth - K1S/hth; print('(5)', float(hth), float(cth)); assert sig(float(cth)) == 1.3e-4
assert float(hc(cth)) - float(hth) < 1e-12
for e in []:
    pass  # 一般式なので check_equation の対象外（\text 付き矢印で表記）
print('THE-ION-076 OK')
