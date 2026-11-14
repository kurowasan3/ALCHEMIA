# 検証スクリプト THE-EQL-062  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-EQL-062/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = "原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸".split()
assert it['field'] in FIELDS, it['field']
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
def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', sorted(eqs))
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    assert p.count('\\begin{enumerate}') <= 1
    n_items = len(re.findall(r'\\item\b', p)) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    assert n_items == len(nums), (n_items, nums)
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from math import log
import numpy as np
from scipy.optimize import brentq
v0 = {40: 2.0, 50: 4.0, 60: 8.0}; tau = {40: None, 50: 10.0, 60: 1.0}
def P(T, t):
    if tau[T] is None: return v0[T]*t
    return v0[T]*tau[T]/log(2)*(1 - 2**(-t/tau[T]))
def Pnum(T, t):
    s = np.linspace(0, t, 200001); y = v0[T]*(2**(-s/tau[T]) if tau[T] else np.ones_like(s))
    return np.trapezoid(y, s)
for T in (50, 60):
    for t in (1, 5, 30): assert abs(P(T,t) - Pnum(T,t)) < 1e-6
print('upper 60', 8.0*1.0/0.693); assert sig(8.0/0.693, 2) == 12
expect = {1: (60, [2.0, 3.9, 5.8]), 5: (50, [10, 17, 11]), 30: (40, [60, 50, 12])}
for t, (Tbest, vals) in expect.items():
    vs = [P(T, t) for T in (40, 50, 60)]
    print(t, [round(v,2) for v in vs])
    assert max((P(T,t), T) for T in (40,50,60))[1] == Tbest
    assert [sig(v, 2) for v in vs] == vals, ([sig(v,2) for v in vs], vals)
# 与えた数値での計算も確認
assert abs(57.7*(1-0.933) - 3.87) < 0.01 and abs(57.7*(1-0.707) - 16.9) < 0.05
x_ap = 0.20/0.693; x_ex = brentq(lambda x: (1-2**(-x))/(x*log(2)) - 0.9, 0.01, 2)
print(x_ap*10, x_ex*10); assert sig(x_ap*10, 2) == 2.9 and sig(x_ex*10, 2) == 3.1 and x_ap < x_ex
print('ALL OK')
