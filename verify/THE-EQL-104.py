# 検証スクリプト（グループ g07）  python3 -I この.py  （単独で動く）
import sys, json, re, os, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

FIELDS = set('''原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡
ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質
有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離
高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'''.split())
FORMATS = {'qa', 'equation', 'compare', 'condition', 'pathway', 'calc', 'structure'}

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

def check_item(path, dmin, dmax):
    it = json.load(open(path, encoding='utf-8'))
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    assert dmin <= it['difficulty'] <= dmax, it['difficulty']
    assert set(it['format'].split(',')) <= FORMATS, it['format']
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    d = os.path.dirname(path)
    for f in ('problem_fig', 'solution_fig'):
        if it[f]: assert os.path.exists(os.path.join(d, it[f])), it[f]
    assert '\\hfill$\\square$' in it['solution']
    return it

def sig(x, n):
    return float(f'{x:.{n}g}')
it = check_item('/home/claude/chemx/todai/outh/items/THE-EQL-104/item.json', 4, 5)
k1 = 0.050; k3 = 0.693/7200; k2 = k1/0.80 - k1 - k3
print('k2', k2, 'k3', k3)
assert sig(k2, 2) == 0.012
xs = lambda f: f*k1/(f*k1 + f*k2 + k3)
print(xs(1), xs(0.5), xs(0.01))
assert sig(xs(0.5), 2) == 0.80 and sig(xs(0.01), 2) == 0.69
def f(x, I=1.0): return I*k1*(1-x) - (I*k2 + k3)*x
def rk4(x, h, I=1.0):
    a = f(x, I); b = f(x + h*a/2, I); c = f(x + h*b/2, I); d = f(x + h*c, I)
    return x + h*(a + 2*b + 2*c + d)/6
x, t, h = 0.0, 0.0, 0.001
target = 0.9*xs(1)
while x < target:
    x = rk4(x, h); t += h
print('t90', t, 'analytic', math.log(10)/(k1+k2+k3), 'with 2.30', 2.30/(k1+k2+k3))
assert sig(t, 2) == 37
x, t, h = 0.80, 0.0, 1.0
while x > 0.10:
    x = rk4(x, h, I=0.0); t += h
print('dark t (h)', t/3600)
assert sig(t/3600, 2) == 6.0
print('THE-EQL-104 OK')
