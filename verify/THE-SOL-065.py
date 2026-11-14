# 検証スクリプト THE-SOL-065  python3 -I この.py
ITEM = '/home/claude/chem6/outh/items/THE-SOL-065/item.json'
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
# ---- THE-SOL-065 の数値検証（SymPy）
import sympy as sp
Kf = sp.Rational(512, 100); M = 122; K = 50
a = sp.symbols('a', positive=True)
def forward(C):
    av = [s for s in sp.solve(sp.Eq(a + 2*K*a**2, C), a) if s > 0][0]
    d = K*av**2; return av, d, av + d
for C, dT_given in ((sp.Rational(6, 100), 0.205), (sp.Rational(2, 10), 0.614)):
    av, d, p = forward(C); dT = Kf*p
    print('C', C, 'a', av, 'd', d, 'dT', float(dT)); assert sig(float(dT), 3) == dT_given
# 逆算（測定値から）
for C, dT, mass, Mexp in ((0.0600, 0.205, 7.32, 183), (0.200, 0.614, 24.4, 203)):
    p = dT/5.12; d = C - p; av = 2*p - C; Kc = d/av**2
    print('Mapp', mass/p, 'K', Kc); assert sig(mass/p, 3) == Mexp and sig(Kc, 2) == 50
assert abs(7.32/122 - 0.0600) < 1e-9 and abs(24.4/122 - 0.200) < 1e-9
av, d, p = forward(sp.Rational(1, 100)); dT = Kf*p
print('C=0.01 a', float(av), 'd', float(d), 'dT', float(dT), 'Mapp', float(sp.Rational(122,100)/p))
assert sig(float(dT), 2) == 0.041
# 極限
Cs = sp.symbols('C', positive=True)
asol = (-1 + sp.sqrt(1 + 8*K*Cs))/(4*K)
Mapp = M*Cs/(asol + K*asol**2)
print('lim0', sp.limit(Mapp, Cs, 0), 'limoo', sp.limit(Mapp, Cs, sp.oo))
assert sp.limit(Mapp, Cs, 0) == 122 and sp.limit(Mapp, Cs, sp.oo) == 244
# 会合度
for C, p in ((0.06, 0.205/5.12), (0.2, 0.614/5.12)):
    print('beta', 2*(C - p)/C)
print('THE-SOL-065 OK')
