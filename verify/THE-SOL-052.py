# 検証スクリプト THE-SOL-052  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-SOL-052/item.json'
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
assert check_all_equations() == 0
r = sp.Rational(585, 945)          # NaCl / NaCl·2H2O
assert abs(float(r) - 0.6190) < 1e-4
# (1)
sol = sp.Rational(200) * sp.Rational(10, 100) / sp.Rational(23, 100); ice = 200 - sol
print('(1)', float(ice), float(sol)); assert sig(float(ice), 3) == 113 and sig(float(sol), 3) == 87.0
# (2)
Dh = 20 / r; ic = 200 - Dh
print('(2)', float(ic), float(Dh)); assert sig(float(Dh), 3) == 32.3 and sig(float(ic), 3) == 168
eut_D = sol * sp.Rational(23, 100) / r; eut_I = sol - eut_D
ratio = eut_I / eut_D; print('共融 氷/二水和物', float(ratio), float(eut_D / sol))
assert sig(float(ratio), 3) == 1.69 and abs(float(eut_D / sol) - 0.372) < 1e-3
# (3)
D = sp.symbols('D')
d1 = sp.solve(sp.Eq(sp.Rational(26), r * D + sp.Rational(248, 1000) * (100 - D)), D)[0]
d0 = sp.solve(sp.Eq(sp.Rational(26), D + sp.Rational(248, 1000) * (100 - D)), D)[0]
print('(3)', float(d1), float(d0)); assert sig(float(d1), 3) == 3.23 and sig(float(d0), 3) == 1.60
assert 26.0 > 24.8 and 26.0 < 26.4   # 20 ℃ では未飽和，-10 ℃ では過飽和
# (4) 場合分け
liq = 100 - sp.Rational(348, 10); nacl_liq = liq * sp.Rational(23, 100)
w_ice = nacl_liq                                  # 氷が先
w_hyd = sp.Rational(348, 10) * r + nacl_liq        # 二水和物が先
print('(4) 氷:', float(w_ice), ' 二水和物:', float(w_hyd))
assert sig(float(w_ice), 3) == 15.0 and w_ice < 23        # 矛盾なし
assert w_hyd > 23 and w_hyd > sp.Rational(264, 10)        # 20 ℃ で溶かせない → 不適
assert sig(float(w_hyd), 3) == 36.5
# 15.0 % の凝固点（表）
assert dict(zip([0, 5, 10, 15, 20, 23], [0, -3.0, -6.6, -10.9, -16.5, -21]))[15] == -10.9
print('THE-SOL-052 OK')
