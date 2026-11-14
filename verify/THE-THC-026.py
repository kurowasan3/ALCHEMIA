# 検証スクリプト THE-THC-026  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-THC-026/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
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
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', b): continue   # 一般式（R・n）は数えない
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n += 1
    print('反応式', n, '本 OK')
    return n
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    assert it['field'] in FIELDS, it['field']
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)

# ---- 個別の検証 ----
import sympy as sp
check_structure()
assert check_all_equations() >= 4
dH = -1273 - (6*(-394) + 6*(-286)); assert dH == 2807
# (2) 係数を連立方程式で
a, b, c, d = sp.symbols('a b c d')
sol = sp.solve([sp.Eq(a, 6), sp.Eq(2*a, 6 + d), sp.Eq(b, 12 + 2*d), sp.Eq(b - c, 0)], [a, b, c, d])
assert sol == {a: 6, b: 24, c: 24, d: 6}
assert check_equation('6CO2 + 24H+ + 24e- -> C6H12O6 + 6H2O') is None
assert check_equation('6CO2 + 12H2O -> C6H12O6 + 6O2 + 6H2O') is None
# (4)(5)
h, cc, NA = sp.Float('6.63e-34'), sp.Float('3.00e8'), sp.Float('6.02e23')
E680 = NA*h*cc/sp.Float('680e-9')/1000
E440 = NA*h*cc/sp.Float('440e-9')/1000
eff680 = 2807/(48*E680)*100; eff440 = 2807/(48*E440)*100
print('E680', E680, 'eff', eff680, ' E440', E440, 'eff', eff440)
assert sig(float(eff680), 3) == 33.2 and sig(float(eff440), 3) == 21.5
# 文献値 ΔfH glucose(s) -1273.3, 光合成の ΔH ≒ +2803〜2808
print('THE-THC-026 OK')
