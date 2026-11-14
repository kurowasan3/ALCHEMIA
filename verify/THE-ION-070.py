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

it = check_item('/home/claude/chem6/out/items/THE-ION-070/item.json')
from sympy import Rational as Q, sqrt, log, nsolve, symbols
kH = Q(32,1000)/Q(10**5)          # mol/(L·Pa)
K1 = Q(5,10**7); Kw = Q(1,10**14)
P = Q(10**5)*Q(4,10000)
C = kH*P
print('C =', float(C)); assert sig(float(C)) == 1.3e-5
h = sqrt(K1*C); pH = float(-log(h,10)); print('pH', pH)
# log2=0.30 を使った値
pH_book = -0.5*(6*0.30-13); print(pH_book); assert abs(pH_book-5.60)<1e-9 and abs(pH-5.60)<0.005
# 水の電離を含めた厳密解（開いた系）： h^2 = K1 C + Kw
assert abs(float(-log(sqrt(K1*C+Kw),10)) - pH) < 1e-3
# 閉じた系とした場合
x = symbols('x', positive=True)
hc = nsolve(x**2 - K1*(C-x), x, 2e-6); print('closed', hc, float(hc/C), float(-log(hc,10)))
assert round(float(-log(hc,10)),2) == 5.64 and 0.15 < float(hc/C) < 0.25
pH2 = float(-log(sqrt(K1*2*C),10)); print(pH2); assert round(pH2,2) == 5.45 and abs(5.60-0.15-5.45)<1e-9
Cs = kH*Q(10**5); hs = sqrt(K1*Cs); pHs = float(-log(hs,10)); print(pHs, float(hs/Cs))
assert round(pHs,2) == 3.90 and abs(-0.5*(4*0.30-9)-3.90)<1e-9
# 第2段階の無視： [CO3^2-]≈K2 ≪ [HCO3-]
K2 = 4.7e-11; assert K2/float(h) < 1e-4
print('THE-ION-070 OK')
