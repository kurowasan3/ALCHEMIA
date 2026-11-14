# 検証スクリプト（グループ g12）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-HYD-033/item.json'
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
a, b = symbols('a b', positive=True)
assert solve(a + b - 2*b, a) == [b]
# 全体の式 = (i)x2 + (ii)x4 + (iii)x1 の原子の確認（X, HCl が消える）
from collections import Counter
def cnt(f, k=1):
    import re
    c = Counter()
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', f.replace('=','')):
        c[el] += k*(int(n) if n else 1)
    return c
L = cnt('C2H4',2)+cnt('Cl2',2)+cnt('C2H4Cl2',4)+cnt('C2H4',2)+cnt('HCl',4)+cnt('O2',1)
Rr = cnt('C2H4Cl2',2)+cnt('C2H3Cl',4)+cnt('HCl',4)+cnt('C2H4Cl2',2)+cnt('H2O',2)
L2 = L - cnt('C2H4Cl2',4) - cnt('HCl',4); R2 = Rr - cnt('C2H4Cl2',4) - cnt('HCl',4)
assert L2 == cnt('C2H4',4)+cnt('Cl2',2)+cnt('O2',1) and R2 == cnt('C2H3Cl',4)+cnt('H2O',2)
Mvc = 2*12 + 3*R('1.0') + R('35.5'); assert Mvc == R('62.5')
n = R(10**6)/Mvc; assert n == 16000
m1 = n/2*71/10**6; m2 = n*71/10**6
print(float(m1), float(m2)); assert sig(float(m1),3) == 0.568 and sig(float(m2),3) == 1.14
assert check_equation('4CH2=CH2 + 2Cl2 + O2 -> 4CH2=CHCl + 2H2O') is None

print('ALL OK')
