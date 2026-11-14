# 検証スクリプト THE-SOL-058  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/THE-SOL-058/item.json'
import sys, json, re
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
    assert it['field'] in FIELDS, it['field']
    assert all(f in ('qa','equation','compare','condition','pathway','calc','structure') for f in it['format'].split(','))
    import os
    for k in ('problem_fig','solution_fig'):
        if it[k]: assert os.path.exists(os.path.join(os.path.dirname(ITEM), it[k])), it[k]
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from sympy import Rational as R, nsolve, symbols, solve, Eq
data = [(R(2), 144), (R(5), 360), (R(7), 504), (R(10), 616), (R(14), 696), (R(18), 776)]
lN, lD = 50, 22
c = symbols('c')
assert all(k == 72*cc for cc, k in data[:3])
assert all(k == 416 + 20*cc for cc, k in data[3:])
cmc = solve(Eq(72*c, 416 + 20*c), c)[0]
assert cmc == 8
assert lN + lD == 72
beta = 1 - R(20, lN)
assert beta == R(3, 5)
# モデル全体で表の値を再現
def kappa(cc):
    if cc <= cmc:
        return (lN + lD)*cc
    return lN*(cc - beta*(cc - cmc)) + lD*cmc
assert all(kappa(cc) == k for cc, k in data)
na18 = 18 - beta*(18 - cmc); mic = (18 - cmc)/64
print('Na', na18, 'micelle', float(mic), 'charge', float(-64*(1 - beta)))
assert na18 == 12 and sig(float(mic), 2) == 0.16 and sig(float(-64*(1-beta)), 2) == -26
# (4)
x = symbols('x', positive=True)
K = cmc**(1 + beta)
print('const', float(K))
assert sig(float(K), 3) == 27.9
res = {}
for cs in (20, 100):
    res[cs] = float(nsolve(x*(x + cs)**beta - K, x, 3))
print(res)
assert sig(res[20], 2) == 4.1 and sig(res[100], 2) == 1.7
assert round(4.0*24.0**0.6, 1) == 26.9 and round(4.2*24.2**0.6, 1) == 28.4
approx = float(K/100**beta)
assert sig(approx, 3) == 1.76 and abs(approx - res[100])/res[100] < 0.015
assert sig(float(K)/101.8**0.6, 3) == 1.74
# 単位換算：1 mmol/L = 1e-6 mol/cm3, λ[S cm2/mol] → κ[S/cm] = λ*1e-6 → μS/cm
assert abs(50*1e-6*1e6 - 50) < 1e-12
import os
s = open(os.path.join(os.path.dirname(ITEM), 'solution.svg'), encoding='utf-8').read()
assert 'viewBox' in s and 'script' not in s

print('ALL OK')
