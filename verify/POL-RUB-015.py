# 検証スクリプト POL-RUB-015  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/POL-RUB-015/item.json'
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
from sympy import symbols, Rational, solve, Eq
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
# 単位の組成：ヒドロキシ酸 − H2O
assert CalcMolFormula(Chem.MolFromSmiles('CC(O)C(=O)O')) == 'C3H6O3'
assert CalcMolFormula(Chem.MolFromSmiles('OCC(=O)O')) == 'C2H4O3'
x = symbols('x')
xs = solve(Eq((2+x)/2, Rational(9,8)), x)[0]; print(xs); assert xs == Rational(1,4)
ok = {k: 1 <= a/b <= 1.5 for k,(a,b) in {'ア':(5,4),'イ':(7,4),'ウ':(4,5),'エ':(11,8)}.items()}
print(ok); assert [k for k,v in ok.items() if not v] == ['イ','ウ']
Mlac, Mgly = 3*12+4*1.0+2*16, 2*12+2*1.0+2*16
assert (Mlac, Mgly) == (72, 58)
Mavg = 0.25*Mlac + 0.75*Mgly; n = 12.3/Mavg
O2 = n*(0.25*3 + 0.75*1.5); print(Mavg, n, O2)
assert abs(O2-0.375) < 1e-9
print(3/72, 1.5/58); assert 3/72 > 1.5/58
assert check_equation('C3H6O3 + 3O2 -> 3CO2 + 3H2O') is None
print('ALL OK')
