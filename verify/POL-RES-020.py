# 検証スクリプト（グループ g10）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/POL-RES-020/item.json'
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
p, q = sp.symbols('p q', positive=True)
c0 = 10
# (1)
r = {}
for K in (4, 400):
    s = [x for x in sp.solve(sp.Eq(p**2, K*(1-p)**2), p) if 0 < x < 1][0]
    r[K] = (s, 1/(1-s)); print('(1)', K, s, 1/(1-s))
assert r[4] == (sp.Rational(2,3), 3) and r[400] == (sp.Rational(20,21), 21)
# (2)
w = {K: K*c0*0.01**2/0.99 for K in (4, 400)}
print('(2)', w, {K: v*18*1000 for K, v in w.items()})
assert sig(w[4], 2) == 0.004 and sig(w[400], 2) == 0.40
assert round(w[4]*18*1000) == 73 and sig(w[400]*18*1000, 2) == 7300
# (3)
qe = [x for x in sp.solve(sp.Eq(4*c0*q**2, 0.040*(1-q)), q) if x > 0][0]
print('(3) exact', qe, 1/qe, 'approx', (4*c0/0.040)**0.5)
assert round(1/qe) == 32 and round((4*c0/0.040)**0.5) == 32
# (4)
qs = [x for x in sp.solve(sp.Eq(4*c0*q**2, (1-q)*(0.040 - c0*(q-0.010))), q) if x > 0][0]
print('(4)', qs, 1/qs, 'water left', 0.040 - c0*(qs-0.01))
assert round(1/qs) == 75 and 0.006 < 0.040 - c0*(qs-0.01) < 0.008
need = c0*(1/32 - 1/100); print('need water', need); assert sig(need, 2) == 0.21
assert '75' in it['answer'] and '32' in it['answer']
print('OK')
