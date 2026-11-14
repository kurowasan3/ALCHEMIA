# 検証スクリプト THE-GAS-047  python3 -I この.py
ITEM = '/home/claude/chem6/outh/items/THE-GAS-047/item.json'
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
# ---- THE-GAS-047 の数値検証（SymPy）
import sympy as sp
x, y, z = sp.symbols('x y z', positive=True)
R = sp.Rational(83, 10)*10**3; T1, T2 = 300, 400; V = sp.Rational(249, 100)
p0 = sp.Rational(1500, 1000); P1 = sp.Rational(1600, 1000)
P2m, P4m, P3m = sp.Rational(1486, 1000), sp.Rational(2333, 1000), sp.Rational(1186, 1000)
pv27 = sp.Rational(36, 1000); pv127 = sp.Rational(246, 100)
a = P1 - p0
n = a*10**5*V/(R*T1)
print('n =', sp.N(n), ' M =', sp.N(sp.Rational(740, 1000)/n))
assert abs(sp.N(n) - 0.0100) < 1e-9 and abs(sp.N(sp.Rational(740,1000)/n) - 74) < 1e-6
# 水がすべて気体なら P2 > p0 のはず
assert P2m < p0
eq1 = sp.Eq(p0 - a*(y/4 - z/2) + pv27, P2m)
eq2 = sp.Eq((p0 + a*(y/4 + z/2))*sp.Rational(T2, T1), P4m)
sol = sp.solve([eq1, eq2], [y, z], dict=True)[0]
print('y,z =', sp.N(sol[y]), sp.N(sol[z]))
assert abs(sol[y] - 6) < sp.Rational(1, 100) and abs(sol[z] - 2) < sp.Rational(1, 100)
yv, zv = 6, 2
w127 = a*yv/2*sp.Rational(T2, T1)
assert w127 < pv127 and a*yv/2 > pv27
eq3 = sp.Eq(p0 - a*(x + sp.Rational(yv, 4) - sp.Rational(zv, 2)) + pv27, P3m)
xv = sp.solve(eq3, x)[0]; print('x =', xv); assert xv == 3
assert 12*3 + 6*1 + 16*2 == 74
# 順方向の再計算で測定値を再現
O2left = p0 - a*(3 + sp.Rational(6, 4) - 1)
assert O2left > 0
assert O2left + a*3 + pv27 == P2m and O2left + pv27 == P3m
print('P4 計算値', sp.N((O2left + a*3 + a*3)*sp.Rational(4, 3)))
assert abs(sp.N((O2left + a*3 + a*3)*sp.Rational(4, 3)) - 2.333) < 1e-3
Vmin = sp.Rational(3, 100)*R*T1/(pv27*10**5)
print('Vmin =', sp.N(Vmin), '->', sig(float(Vmin), 2)); assert sig(float(Vmin), 2) == 21
print('THE-GAS-047 OK')
