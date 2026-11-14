# 検証スクリプト ORG-CBX-045（グループ g10）  python3 -I この.py
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

FIELDS = set("原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸".split())

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
    return it

it = check_item('/home/claude/chem6/outh/items/ORG-CBX-045/item.json')
from collections import Counter
from fractions import Fraction as F
# (2) 機構 I のシミュレーション（段階を少しずつ適用）
s = Counter({'C2H5OH': F(100), 'CH3I': F(10), 'H2O': F(20)})  # 単位 0.01 mol
step = F(1, 10)
def run(s):
    for _ in range(100000):
        moved = False
        for R, alc in (('CH3', 'CH3OH'), ('C2H5', 'C2H5OH')):
            # ① ROH + HI -> RI + H2O
            x = min(step, s[alc], s['HI'])
            if x > 0: s[alc] -= x; s['HI'] -= x; s[R+'I'] += x; s['H2O'] += x; moved = True
            # ② RI + CO -> RCOI
            x = min(step, s[R+'I'])
            if x > 0: s[R+'I'] -= x; s[R+'COI'] += x; moved = True
            # ③ RCOI + H2O -> RCOOH + HI
            x = min(step, s[R+'COI'], s['H2O'])
            if x > 0: s[R+'COI'] -= x; s['H2O'] -= x; s[R+'COOH'] += x; s['HI'] += x; moved = True
        if not moved: return s
run(s)
res = {k: float(v/100) for k, v in s.items() if v}
print('機構 I 終了時', res)
assert res == {'CH3COOH': 0.1, 'C2H5COOH': 1.0, 'HI': 0.1, 'H2O': 0.1}
# (3) 原子収支: CH3COOCH3 1.00 + CO 1.00 + H2O 0.050 -> a (CH3CO)2O + b CH3COOH
from sympy import symbols, solve, Rational as Q
a, b = symbols('a b')
sol = solve([4*a + 2*b - (3*1 + 1), 3*a + 2*b - (2*1 + 1 + Q(5,100)), 6*a + 4*b - (6*1 + 2*Q(5,100))], [a, b], dict=True)
print(sol); assert sol == [{a: Q(95,100), b: Q(1,10)}]
for e in ['C2H5OH + HI -> C2H5I + H2O','C2H5I + CO -> C2H5COI','C2H5COI + H2O -> C2H5COOH + HI','CH3COOCH3 + CO + H2O -> 2CH3COOH','C2H5OH + CO -> C2H5COOH']:
    assert check_equation(e) is None, e
print('OK')
