# 検証スクリプト THE-SOL-061  python3 -I この.py
ITEM = '/home/claude/chem6/out/items/THE-SOL-061/item.json'
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
# ---- THE-SOL-061 の数値検証（SymPy）
import sympy as sp
ns = sp.Rational(684, 100)/342; assert ns == sp.Rational(2, 100)
w = sp.Rational(667, 10)
n = ns*(200 + w)/(200 - w); print('n =', sp.N(n)); assert sig(float(n), 3) == 0.0400
M = sp.Rational(117, 100)/n; print('M_app =', sp.N(M)); assert sig(float(M), 2) == 29
print('式量/見かけ =', sp.N(sp.Rational(585, 10)/M)); assert abs(float(sp.Rational(585,10)/M) - 2) < 0.02
# 正しい NaCl 0.0200 mol (粒子0.0400) なら移動量は？
W = sp.symbols('W'); wtrue = sp.solve(sp.Eq(ns/(200 - W), sp.Rational(4,100)/(200 + W)), W)[0]
print('NaCl 完全電離での移動量', sp.N(wtrue)); assert abs(float(wtrue) - 66.7) < 0.05
n4 = ns*240/160; M4 = sp.Rational(117,100)/n4; print('打ち切り M =', sp.N(M4)); assert sig(float(M4), 2) == 39 and M4 > M
ng = sp.Rational(117,100)/180; assert ng < ns
wg = sp.solve(sp.Eq(ns/(200 + W), ng/(200 - W)), W)[0]; print('グルコース移動量', sp.N(wg)); assert sig(float(wg), 2) == 100
# 等モル分率（厳密なラウール）との比較：近似の妥当性
x = sp.symbols('x')
wx = sp.nsolve(ns/(ns + (200 - x)/18) - sp.Rational(4,100)/(sp.Rational(4,100) + (200 + x)/18), x, 66)
print('モル分率で解いた移動量', wx); assert abs(wx - wtrue) < 0.5
print('THE-SOL-061 OK')
