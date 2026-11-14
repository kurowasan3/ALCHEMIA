# 検証スクリプト（グループ g06）  python3 -I この.py
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
    for f in it['format'].split(','):
        assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    assert it['answer'].count('$') % 2 == 0 and it['solution'].count('$') % 2 == 0 and it['problem'].count('$') % 2 == 0
    return it

def sig(x, n=2):
    from math import log10, floor
    if x == 0: return 0
    d = n - 1 - floor(log10(abs(x)))
    return round(x, d)

it = check_item('/home/claude/chem6/outh/items/THE-ELC-050/item.json')
F = 9.65e4
Ed_Cu = 1.23-0.34; Ed_H = 1.23-0.00
assert round(Ed_Cu,2) == 0.89 and Ed_H == 1.23
assert 1.0 > Ed_Cu and 1.0 < Ed_H
I0 = (2.0-Ed_Cu)/0.37; t = 0.030*2*F/I0; I1 = (2.0-Ed_H)/0.37
print(I0, t, I1); assert sig(I0) == 3.0 and sig(t) == 1.9e3 and sig(I1) == 2.1
I4 = (6.0-Ed_H)/10; t4 = 0.10*F/I4; t4n = 0.10*F/0.60
print(I4, t4, t4n, 1-t4n/t4, Ed_H/6.0)
assert sig(t4) == 2.0e4 and abs(100*(1-t4n/t4) - 20.5) < 1e-9 and abs(100*Ed_H/6.0 - 20.5) < 1e-9
assert abs((1-t4n/t4) - Ed_H/6.0) < 1e-12
for e in ['2Cu + O2 + 4H+ -> 2Cu^2+ + 2H2O','2H2O -> O2 + 4H+ + 4e-','Cu -> Cu^2+ + 2e-','O2 + 4H+ + 4e- -> 2H2O']:
    assert check_equation(e) is None
print('THE-ELC-050 OK')
