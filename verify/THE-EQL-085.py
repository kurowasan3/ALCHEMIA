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

it = check_item('/home/claude/chem6/outh/items/THE-EQL-085/item.json')
assert it['difficulty'] in (4, 5)
from sympy import symbols, nsolve, Rational as Q, solve, cbrt, sqrt, N
y = symbols('y', positive=True)
Kp = Q(4,10**13); p0 = Q(10**5)
# (1) 厳密：(2y)^2 y / (p0-2y)^2 = Kp
y1 = nsolve((2*y)**2*y - Kp*(p0-2*y)**2, y, 0.1); print('(1) pH2', 2*y1); assert sig(float(2*y1)) == 0.20
# (2) 厳密：p_O2 = 1000 + y
y2 = nsolve((2*y)**2*(1000+y) - Kp*(p0-2*y)**2, y, 1e-3); print('(2) pH2', 2*y2); assert sig(float(2*y2)) == 2.0e-3
# (3) 容積2倍
y1b = nsolve((2*y)**2*y - Kp*(p0/2-2*y)**2, y, 0.06)
y2b = nsolve((2*y)**2*(500+y) - Kp*(p0/2-2*y)**2, y, 1e-3)
r1 = (2*y1b/(p0/2))/(2*y1/p0); r2 = (2*y2b/(p0/2))/(2*y2/p0)
print('(3)', r1, r2); assert round(float(r1),2) == 1.26 and round(float(r2),2) == 1.41
# (4)
y4 = nsolve((2*y)**2*(Q(1,10)+y) - Kp*(p0-2*y)**2, y, 0.075); print('(4) y', y4, 2*y4); assert sig(float(2*y4)) == 0.15
f = lambda v: v**3 + 0.10*v**2 - 1.0e-3
print(f(0.075), f(0.076)); assert f(0.075) < 0 < f(0.076)
approx_b = (4.0e-3/0.10)**0.5; print('approx', approx_b); assert abs(approx_b-0.20) < 1e-12
print('ずれ', 0.20/float(2*y4)-1)
# (5)
for pa, expect in [(1e-3, 0.199), (1e-2, 0.19), (10, 0.020), (1e3, 0.0020)]:
    v = nsolve((2*y)**2*(pa+y) - 4.0e-3, y, 0.1 if pa < 1 else (1e-3/pa)**0.5)
    print(pa, 2*v); assert abs(float(2*v)-expect)/expect < 0.03
# 約1150 K の妥当性：ΔG°(2H2O→2H2+O2) ≈ 483.6 − 0.0888 T kJ, Kp(bar)=exp(−ΔG/RT), Pa に換算
import math
T = 1150; K_bar = math.exp(-(483600-88.8*T)/(8.314*T)); print('Kp(Pa) at 1150K ≈', K_bar*1e5)
assert 1e-13 < K_bar*1e5 < 2e-12
assert check_equation('2H2O <=> 2H2 + O2') is None
print('THE-EQL-085 OK')
