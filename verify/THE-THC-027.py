# 検証スクリプト THE-THC-027  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-THC-027/item.json'
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
check_all_equations()
x, y, n = sp.symbols('x y n')
s = sp.solve([6*x + y - (2*717 + 3*436 + 84), 8*x + 2*y - (3*717 + 4*436 + 105)], [x, y])
assert s == {x: 413, y: 348}, s
f = sp.expand(717*n + 436*(n + 1) - (413*(2*n + 2) + 348*(n - 1)))
assert f == -21*n - 42
assert sp.solve(sp.Eq(f, -168), n) == [6]
assert f.subs(n, 2) == -84 and f.subs(n, 3) == -105 and f.subs(n, 4) == -126
c = sp.expand(-394*n - 286*(n + 1) - f); assert c == -659*n - 244
# O2 係数の確認（n=1..6）：原子のつり合い
for k in range(1, 7):
    o2 = sp.Rational(3*k + 1, 2)
    assert 2*o2 == 2*k + (k + 1)
# 異性体の結合数（RDKit）
from rdkit import Chem
for smi in ('CCCC', 'CC(C)C'):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ch = sum(1 for b in m.GetBonds() if {b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()} == {'C', 'H'})
    cc = sum(1 for b in m.GetBonds() if b.GetBeginAtom().GetSymbol() == b.GetEndAtom().GetSymbol() == 'C')
    assert (ch, cc) == (10, 3)
assert -126 - (-135) == 9
# 文献値（kJ/mol）：C2H6 -84.0, C3H8 -104.7, butane -125.6, isobutane -134.2, hexane(g) -166.9
print('THE-THC-027 OK')
