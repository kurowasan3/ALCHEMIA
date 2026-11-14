# 検証スクリプト THE-EQL-057  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-EQL-057/item.json'
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
import math
import sympy as sp
check_structure()
check_all_equations()
R0, P0, TS, TSc = 0, 50, 180, 120
assert (P0 - R0, TS - R0, TS - P0) == (50, 180, 130)
assert (TSc - R0, TSc - P0) == (120, 70)
assert (TS - R0) == (TS - P0) + (P0 - R0)
# (4) マクスウェル–ボルツマン分布（3次元，エネルギー）で T2>T1 のとき：面積一定・山が低く右へ・E2 以上の割合が増える
def f(E, T): return 2 * math.sqrt(E / math.pi) * T ** -1.5 * math.exp(-E / T)
def integ(g, a, b, n=200000):
    h = (b - a) / n
    return sum(g(a + (i + 0.5) * h) for i in range(n)) * h
T1, T2 = 1.0, 1.3
A1, A2 = integ(lambda E: f(E, T1), 0, 60), integ(lambda E: f(E, T2), 0, 60)
assert abs(A1 - 1) < 1e-3 and abs(A2 - 1) < 1e-3
peak1 = (T1 / 2, f(T1 / 2, T1)); peak2 = (T2 / 2, f(T2 / 2, T2))
assert peak2[0] > peak1[0] and peak2[1] < peak1[1]
E2 = 3.6
assert integ(lambda E: f(E, T2), E2, 60) > integ(lambda E: f(E, T1), E2, 60)
# (5)
e = sp.Float(180e3) / sp.Float(6.02e23)
print('(5)', e); assert sig(float(e), 2) == 3.0e-19
print('THE-EQL-057 OK')
