# 検証スクリプト POL-FIB-013  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/POL-FIB-013/item.json'
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
from sympy import symbols, Eq, solve, Rational
x, n = symbols('x n', positive=True)
H, C, O = 1.0, 12, 16
u = 2*C + 4*H + O                      # C2H4O = 44
acetal_pair = 5*C + 8*H + 2*O          # C5H8O2 = 100
assert u == 44 and acetal_pair == 100
assert check_equation('2C2H4O + CH2O -> C5H8O2 + H2O') is None
sol = solve([Eq(n*(44+6*x), Rational(928,10)), Eq(n*(2+x/2), Rational(1936,10)/44)], [x, n], dict=True)
print(sol); s = sol[0]
assert abs(float(s[x]) - 0.40) < 1e-9 and abs(float(s[n]) - 2.00) < 1e-9
hcho = float(s[n]*s[x]/2); m = hcho*30
print('HCHO mol', hcho, 'g', m, 'formalin', m/0.30)
assert abs(m/0.30 - 40) < 1e-9 and sig(m/0.30, 2) == 40
fC = lambda xx: 12*(2+xx/2)/(44+6*xx)
print(fC(0), fC(0.4), fC(1))
assert fC(0) < fC(0.4) < fC(1)
# 単調増加（導関数の符号）
import sympy as sp
X = sp.symbols('X'); d = sp.diff(12*(2+X/2)/(44+6*X), X); print(sp.simplify(d))
assert all(float(d.subs(X, v)) > 0 for v in [0, 0.3, 0.7, 1])
print('ALL OK')
