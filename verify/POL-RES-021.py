# 検証スクリプト（グループ g12）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/POL-RES-021/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
assert it['field'] in FIELDS, it['field']

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import symbols, solve, Rational as R
E, V, A, Vac = 28, 86, 44, 4*12+6+2*16
assert Vac == 86 and 2*12+4+16 == A and 2*12+4 == E and V - A == 42
a, v = symbols('a v')
sol = solve([42*a - (R('9.90') - R('8.64')), v - a - (R('0.100') - R('1.00')*R('0.0800'))], [a, v])
print(sol); assert sol[a] == R('0.03') and sol[v] == R('0.05')
mE = R('9.90') - sol[v]*V; nE = mE/E; assert nE == R('0.2')
frac = sol[v]/(sol[v]+nE); assert frac == R('0.2')
deg = sol[a]/sol[v]*100; assert deg == 60
assert nE*E + sol[v]*A == R('8.64') - (sol[v]-sol[a])*42 == R('7.80')
# (5) メチルオレンジ：終点 pH 3.1〜4.4 での酢酸イオンのプロトン化率（pKa 4.76）
import math
for pH in (3.1, 4.0, 4.4):
    f = 1/(1+10**(pH-4.76)); vol = 80.0 + 20.0*f; print('pH', pH, 'f', round(f,3), 'V', round(vol,1))
    assert 93.5 < vol <= 100
# 誤った結論の数値
nE_wrong = (R('9.90') - R('0.03')*86)/28; print(float(nE_wrong), float(R('0.03')/(R('0.03')+nE_wrong)))
assert sig(float(nE_wrong),3) == 0.261 and abs(float(R('0.03')/(R('0.03')+nE_wrong)) - 0.10) < 0.005
for e in ['CH3COOCH(CH3)2 + CH3OH -> CH3COOCH3 + (CH3)2CHOH','CH3COOCH(CH3)2 + NaOH -> CH3COONa + (CH3)2CHOH','CH3COO- + H+ -> CH3COOH']:
    assert check_equation(e) is None

print('ALL OK')
