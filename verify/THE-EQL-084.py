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
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    print('小問数', np_, '答え', na)
    return it

it = check_item('/home/claude/chem6/outh/items/THE-EQL-084/item.json')
import math
from sympy import symbols, solve, Rational, simplify, Eq
h, Ka = symbols('h K_a', positive=True)
k1, k2 = Rational(2,10**5), Rational(6,10**4)
kobs = (k1*h + k2*Ka)/(h + Ka)
Kas = solve(Eq(kobs.subs(h, Rational(1,10**4)), Rational(31,10**5)), Ka); print('Ka', Kas); assert Kas == [Rational(1,10**4)]
kk = kobs.subs(Ka, Rational(1,10**4))
hs = solve(Eq(kk, Rational(1,10**4)), h); print('h', hs, -math.log10(float(hs[0]))); assert hs == [Rational(625,10**6)]
assert round(2 + 4*0.30, 1) == 3.2 and round(-math.log10(float(hs[0])), 1) == 3.2
h3 = Rational(1,1000); frac = Rational(1,10**4)/(h3+Rational(1,10**4)); share = k2*Rational(1,10**4)/(k1*h3 + k2*Rational(1,10**4))
print(float(frac), float(share)); assert f'{float(frac)*100:.2g}' == '9.1' and float(share) == 0.75
# 機構の区別ができないこと
X, kp, Kw = symbols('X kp K_w', positive=True)
HX = h*X/Ka; OH = Kw/h
assert simplify(kp*HX*OH - kp*Kw/Ka*X) == 0
kpv = k2*Rational(1,10**4)/Rational(1,10**14); print("k'", kpv); assert kpv == 6*10**6
for pH in range(1, 8):
    hh = 10.0**-pH; ka = 1e-4; c = 1.0
    HXc = c*hh/(hh+ka); Xc = c*ka/(hh+ka)
    va = 2e-5*HXc + 6e-4*Xc; vb = 2e-5*HXc + 6e6*HXc*(1e-14/hh)
    assert math.isclose(va, vb, rel_tol=1e-12)
print('THE-EQL-084 OK')
