ITEM = '/home/claude/chemx/x1/outh/items/INO-PPT-052/item.json'
# 検証スクリプト（x1 g05）  python3 -I この.py
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
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS, it['field']
for f in it['format'].split(','):
    assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
import os
d = os.path.dirname(ITEM)
for k in ('problem_fig','solution_fig'):
    if it[k]: assert os.path.exists(os.path.join(d, it[k])), it[k]

# ---- 個別の検証 ----
from sympy import symbols, diff, solve, sqrt, log
import numpy as np
K1, K2, Kb, Kw, N = 2.0e3, 8.0e3, 1.8e-5, 1.0e-14, 0.100
L = symbols('L', positive=True)
f1 = K1*L/(1+K1*L+K1*K2*L**2)
Lm = solve(diff(f1, L), L)[0]
chk('Lmax', float(Lm), 2.5e-4, 2); chk('f1max', float(f1.subs(L, Lm)), 0.20, 2)
L2 = 1/K2; chk('L(A1=A2)', L2, 1.3e-4, 2)
H = Kw/Kb*(N-L2)/L2; chk('H', H, 4.4e-7, 2); chk('pH', -math.log10(H), 6.4, 2)
def a1_ever_max(k1, k2):
    for x in np.logspace(-9, 1, 4000):
        a0, a1, a2 = 1, k1*x, k1*k2*x*x
        if a1 > a0 and a1 > a2: return True
    return False
for k1, k2 in [(2e3, 8e3), (8e3, 2e3), (1e3, 1e3*0.99), (1e3, 1.01e3), (5e3, 5e3)]:
    assert a1_ever_max(k1, k2) == (k1 > k2), (k1, k2)
print('condition K1>K2 confirmed')

print('ALL OK')
