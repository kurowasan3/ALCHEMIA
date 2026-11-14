# 検証スクリプト THE-SOL-060  python3 -I この.py
ITEM = '/home/claude/chem6/out/items/THE-SOL-060/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = set('原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split())
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
def check_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                r = check_equation(b)
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    return n
def check_structure():
    assert it['field'] in FIELDS, it['field']
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
check_equations()
# ---- THE-SOL-060 の数値検証（SymPy）
import sympy as sp
x, y = sp.symbols('x y')
fC, fW = sp.Rational(160, 250), sp.Rational(90, 250)
assert 160 + 5*18 == 250
SC, SB = sp.Rational(20, 100), sp.Rational(5, 100)
C0, B0, W0 = 45, 6, 150
y1 = sp.solve(sp.Eq(C0 - fC*y, SC*(W0 - fW*y)), y)[0]; print('y1', float(y1)); assert sig(float(y1), 3) == 26.4
assert W0 - fW*y1 > B0/SB and B0/SB == 120
sol = sp.solve([sp.Eq(W0 - x - fW*y, B0/SB), sp.Eq(C0 - fC*y, SC*(W0 - x - fW*y))], [x, y]); print({k: float(v) for k, v in sol.items()})
assert sig(float(sol[x]), 3) == 18.2 and sig(float(sol[y]), 3) == 32.8
y3 = sp.solve(sp.Eq(C0 - fC*y, SC*(W0 - 30 - fW*y)), y)[0]; W3 = W0 - 30 - fW*y3
dep = B0 - SB*W3; print('y3', float(y3), 'W3', float(W3), 'H3BO3析出', float(dep)); assert sig(float(dep), 2) == 0.67
# 60℃で溶けていること（与えた条件の確認：60℃の CuSO4 溶解度は約 40 g/100 g水、H3BO3 約 15 g）
assert 45/150 < 0.40 and 6/150 < 0.15
print('THE-SOL-060 OK')
