# 検証スクリプト THE-EQL-086（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

FIELDS = set('原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split())

def check_item(path, dmin, dmax):
    it = json.load(open(path, encoding='utf-8'))
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
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    assert dmin <= it['difficulty'] <= dmax
    for f in it['format'].split(','):
        assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    nums = [int(m[1]) for m in re.findall(r'(^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, np_ + 1)), ('小問の数', np_, nums)
    assert '\\hfill$\\square$' in it['solution']
    return it

it = check_item('/home/claude/chem6/outh/items/THE-EQL-086/item.json', 4, 5)

from sympy import symbols, solve, Eq
x = (1.0e5-2.4e4)/2.4e4; print(x); assert round(x,1) == 3.2
n = 2.4e4*10.0/(8.31e3*1073); print(n); assert round(n,3) == 0.027
m = symbols('m')
m800 = solve(Eq(m+1, 0.24*(11*m+1)), m)[0]; m850 = solve(Eq(m+1, 0.50*(11*m+1)), m)[0]
print(m800, m850); assert round(float(m800),2) == 0.46 and round(float(m850),2) == 0.11
mh = 178/802; print(mh); assert round(mh,2) == 0.22
assert m800 > mh and m850 < mh
xf = (mh+1)/(11*mh+1); print(xf); assert xf < 0.50

print('ALL OK')
