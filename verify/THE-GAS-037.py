# 検証スクリプト THE-GAS-037  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-GAS-037/item.json'
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
R = sp.Rational(831, 100) * 10**3
vp = {90: 0.701e5, 95: 0.845e5, 100: 1.013e5, 110: 1.433e5, 120: 1.985e5, 130: 2.701e5}
# (1) 沸点 = 飽和蒸気圧が大気圧に等しい温度
P0 = 8.45e4
t1 = [t for t, p in vp.items() if abs(p - P0) < 1]
assert t1 == [95]
# (2) 気泡中の水蒸気
n = sp.Float(P0) * sp.Float(0.60e-3) / (R * 368)
print('(2) n =', n); assert sig(float(n), 2) == 1.7e-05
# (3) 1.433e5 Pa -> 110 ℃
assert [t for t, p in vp.items() if abs(p - 1.433e5) < 1] == [110]
# (5) 空気の分圧（定積）＋飽和蒸気圧
p_air = sp.Float(8.45e4) * sp.Rational(383, 293)
P = p_air + sp.Float(1.433e5)
print('(5) p_air =', p_air, ' P =', P); assert sig(float(P), 2) == 2.5e5
# 実測の水の飽和蒸気圧との整合（文献値 kPa：90→70.1, 95→84.5, 100→101.3, 110→143.3, 120→198.5, 130→270.1）
lit = {90: 70.1, 95: 84.5, 100: 101.3, 110: 143.3, 120: 198.5, 130: 270.1}
for t in vp: assert abs(vp[t] / 1e3 - lit[t]) < 0.2
print('THE-GAS-037 OK')
