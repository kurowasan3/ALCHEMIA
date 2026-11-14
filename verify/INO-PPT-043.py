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
it = check_item('/home/claude/chemx/todai/outh/items/INO-PPT-043/item.json', 4, 5)
KK, KN = 1.0e6, 1.0e4
def bisect(f, lo, hi, n=200):
    flo = f(lo)
    for _ in range(n):
        mid = (lo + hi)/2; fm = f(mid)
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return (lo + hi)/2
def solve(At, Kt, Nt):
    # unknown free A = c ; bound K = KK*c*Kt/(1+KK*c) ; bound Na = KN*c*Nt/(1+KN*c)
    g = lambda c: c + KK*c*Kt/(1 + KK*c) + KN*c*Nt/(1 + KN*c) - At
    c = bisect(g, 0.0, At)
    yK = KK*c*Kt/(1 + KK*c); zN = KN*c*Nt/(1 + KN*c)
    return c, yK, zN
c, y, z = solve(1e-3, 1e-3, 0.0); f1 = y/1e-3; print('(1)', f1)
assert sig(f1*100, 2) == 97
c, y, z = solve(1e-3, 1e-3, 0.100); f2 = y/1e-3; print('(2)', f2, 'ANa', z, z/0.1)
assert sig(f2*100, 2) == 38 and z/0.1 < 0.01
print('(3)', y/z); assert sig(y/z, 2) == 0.62
# approximate route from solution text
u = (-1 + math.sqrt(1 + 4*999*1e-3))/1998; print('approx u', u, math.sqrt(4.996), math.sqrt(4.004e-9))
assert sig(u, 3) == 6.18e-4
At = bisect(lambda A: solve(A, 1e-3, 0.100)[1]/1e-3 - 0.90, 1e-3, 1e-1)
print('(4)', At)
assert sig(At, 2) == 9.2e-3
z_approx = KN*0.100*9.0e-6; print('approx total', 9e-6 + 9e-4 + z_approx)
assert sig(9e-6 + 9e-4 + z_approx, 2) == 9.9e-3
print('INO-PPT-043 OK')
