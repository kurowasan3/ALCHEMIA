# 検証スクリプト THE-SOL-050  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-SOL-050/item.json'
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
x = sp.symbols('x', positive=True)
x1 = sp.solve(sp.Eq(x / (100 + x), sp.Rational(15, 100)), x)[0]
print('(1)', float(x1)); assert sig(float(x1), 2) == 18
melt = sp.Rational(50, 10) / sp.Rational(15, 100) - 5
print('(2)', float(melt)); assert sig(float(melt), 2) == 28
w5 = sp.Rational(30, 130) * 100
# 図の飽和曲線（折れ線 (23,-21),(24.8,-10),(26.3,0)）で -10.9 ℃ の飽和濃度を内挿
def sat(t):
    pts = [(-21, 23), (-10, 24.8), (0, 26.3)]
    for (t1, w1), (t2, w2) in zip(pts, pts[1:]):
        if t1 <= t <= t2: return w1 + (w2 - w1) * (t - t1) / (t2 - t1)
s = sat(-10.9)
print('(5) w =', float(w5), ' 飽和', s); assert 15 < float(w5) < s
assert x1 < 30   # 氷はすべて融ける
# 凝固点の表の単調性（濃いほど凝固点が低い）
fp = [0, -3.0, -6.6, -10.9, -16.5, -21]
assert all(a > b for a, b in zip(fp, fp[1:]))
print('THE-SOL-050 OK')
