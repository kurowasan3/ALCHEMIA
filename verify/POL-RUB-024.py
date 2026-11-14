# 検証スクリプト（グループ g10）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/POL-RUB-024/item.json'
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
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS, it['field']
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert it['solution'].rstrip().endswith('\\hfill$\\square$')

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import sympy as sp
l = sp.symbols('l', positive=True)
g = lambda x: x - 1/x**2
T = 300 * g(2.00) / g(1.80); print('(1)', T); assert sig(T, 3) == 352
alpha = 2.2e-4
lstar = [x for x in sp.solve(sp.Eq((l**3-1)/(l**3+2), alpha*300), l) if x.is_real][0]; print('(2)', lstar)
assert sig(float(lstar), 3) == 1.07
def Lratio(l1, T1=300, T2=301):
    t = T1 * g(l1) / T2
    l2 = sp.nsolve(l - 1/l**2 - t, l, l1)
    return float(l2 / l1 * (1 + alpha*(T2 - T1)))
for l1, up in [(1.05, True), (1.06, True), (1.07, False), (1.10, False), (2.0, False)]:
    r = Lratio(l1); print(l1, r); assert (r > 1) == up
Wint = sp.integrate(l - 1/l**2, (l, 1, 3)); assert Wint == sp.Rational(10, 3)
W = float(Wint) * 4.0e-4 * 8.31 * 300; print('(3)', W); assert sig(W, 2) == 3.3
dT = W / 4.0; print('(4)', dT, -W/300); assert sig(dT, 2) == 0.83 and sig(-W/300, 2) == -0.011
for s in ['352', '1.07', '3.3', '0.83']:
    assert s in it['answer'], s
print('OK')
