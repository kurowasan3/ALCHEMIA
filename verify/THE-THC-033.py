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

def check_item(path):
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
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    print('小問数', np_, '答え', na)
    return it

it = check_item('/home/claude/chem6/out/items/THE-THC-033/item.json')
from sympy import symbols, solve, Rational, Eq
for e in ['C + H2O -> CO + H2', 'CH4 + H2O -> CO + 3H2', 'CO + 2H2 -> CH3OH', 'CO + H2O -> CO2 + H2']:
    assert check_equation(e) is None, e
# 生成エンタルピーからの確認
f = dict(CH4=-74.8, H2O=-241.8, CO=-110.5, CO2=-393.5, CH3OH=-201.0)
print('a', f['CO']-f['H2O'], 'b', f['CO']-f['CH4']-f['H2O'], 'c', f['CH3OH']-f['CO'], 'd', f['CO2']-f['CO']-f['H2O'])
# (1)
print('(a) CO 余り', 1 - Rational(1,2), '(b) H2 余り', 3 - 2)
# (2)
m = 1.2e3/12*Rational(2,3)*32; print('(2) g', float(m)); assert round(float(m)/1000, 1) == 2.1
# (3)
x, y = symbols('x y'); s = solve(Eq(x+3*y, 2*(x+y)), x); print('(3)', s); assert s == [y]
# (4)
H2 = 3*131 - 41 + 2*(-91); H3 = 131 + 206 + 2*(-91); print('(4)', H2, H3); assert (H2, H3) == (170, 155)
print('THE-THC-033 OK')
