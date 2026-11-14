# 検証スクリプト THE-ELC-052（グループ g07）  python3 -I この.py
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

it = check_item('/home/claude/chem6/outh/items/THE-ELC-052/item.json', 4, 5)

from sympy import symbols, diff, solve, log, Matrix, Rational
f = symbols('f'); assert Rational(2,3) in solve(diff(f**2*(1-f), f), f)
pH = -log(2*3.0e-8, 10).evalf(); print('pH', pH); assert round(float(pH),1) == 7.2
# 経路Ⅱ：ClO3- 4 個あたり Cl2 12 個 (電子 24) + 陽極 12 電子
assert (24+12)/4 == 9
# 経路Ⅰ・経路Ⅱ・水の酸化の (電子, ClO3-, O2, H+) の係数（ClO3- 1 mol 又は O2 1 mol あたり）
# 経路Ⅰ: e 6, ClO3 1, O2 0, H+ 0（全体 NaCl+3H2O→NaClO3+3H2 で H+ 正味0）
# 経路Ⅱ: e 9, ClO3 1, O2 3/4, H+ 3（経路Ⅰより余分）
# 水の酸化: e 4, ClO3 0, O2 1, H+ 4
M = Matrix([[6, 9, 4], [1, 1, 0], [0, Rational(3,4), 1], [0, 3, 4]])
print('rank', M.rank()); assert M.rank() == 2   # 内訳は決まらない
n = (9.00 - 4*0.150)/6; print(n, 6*n/9.00); assert abs(n-1.40)<1e-9 and round(6*n/9.00*100,1) == 93.3

print('ALL OK')
