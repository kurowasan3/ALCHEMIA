# 検証スクリプト（京大 g04）  python3 -I この.py
import sys, json, re, math
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

FIELDS = ['原子の構造と電子配置','元素の周期表と周期性','化学結合と結晶','物質量と化学反応式','酸と塩基','酸化還元反応','電池と電気分解','気体の性質','溶液の性質','化学反応とエネルギー','反応速度と化学平衡','電離平衡と溶解平衡','ハロゲン・貴ガス・水素','酸素・硫黄','窒素・リン','炭素・ケイ素','アルカリ金属・アルカリ土類金属','アルミニウム・亜鉛・スズ・鉛','遷移元素（鉄・銅・銀・クロム・マンガン）','金属イオンの沈殿・錯イオン・系統分析','気体の製法と性質','有機化合物の構造と元素分析','脂肪族炭化水素','アルコール・アルデヒド・ケトン','カルボン酸・エステル・油脂','芳香族化合物','有機化合物の検出反応と分離','高分子の基礎と合成樹脂','合成繊維','ゴム・機能性高分子','糖類','アミノ酸・タンパク質・核酸']

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
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it

def sig(x, n=2):
    return float(f'{x:.{n-1}e}')

it = check_item('/home/claude/chemx/kyodai/out/items/POL-SUG-043/item.json')
free = 0.100*16.0e-3
# 収支をたどる
x = 1.20e-3   # 仮のエステル量
naoh = 0.500*20.0e-3; hcl = 0.500*20.0e-3
excess = naoh - x; left_hcl = hcl - excess   # 過剰NaOHを中和した後に余る HCl
assert abs(left_hcl - x) < 1e-15
ester = 0.100*12.0e-3
DE = ester/(free+ester)*100; assert round(DE)==43
Ka = 10**-3.5
f35 = free*Ka/(Ka+10**-3.5)/(free+ester)*100; f55 = free*Ka/(Ka+10**-5.5)/(free+ester)*100
assert round(f35)==29 and round(f55)==57, (f35,f55)
M={'C':12.011,'H':1.008,'O':15.999}
fm=lambda c,h,o: c*M['C']+h*M['H']+o*M['O']
print('C6H8O6', fm(6,8,6), 'C7H10O6', fm(7,10,6))
assert round(fm(6,8,6))==176 and round(fm(7,10,6))==190
mass = free*176 + ester*190; assert round(mass/0.800*100)==64
print(free, ester, DE, f35, f55, mass)
print('OK')
