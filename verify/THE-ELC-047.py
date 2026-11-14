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

it = check_item('/home/claude/chem6/out/items/THE-ELC-047/item.json')
E = {'Mg': -2.37, 'Al': -1.66, 'Zn': -0.76, 'Fe2/Fe': -0.44, 'H': 0.0, 'Cu': 0.34, 'I2': 0.54, 'Fe3/Fe2': 0.77, 'Ag': 0.80, 'Br2': 1.09, 'Cl2': 1.36}
assert round(E['Ag']-E['Mg'],2) == 3.17 and round(E['Cu']-E['Al'],2) == 2.00 and round(E['Ag']-E['Cu'],2) == 0.46
# (2) 酸化剤側, 還元剤側
rx = {'ア': ('Fe3/Fe2','I2'), 'イ': ('Fe3/Fe2','Br2'), 'ウ': ('Fe3/Fe2','Cu'), 'エ': ('Cu','Ag'), 'オ': ('Fe3/Fe2','Fe2/Fe'), 'カ': ('Br2','Cl2')}
go = [k for k,(ox,red) in rx.items() if E[ox] > E[red]]
print(go); assert go == ['ア','ウ','オ']
ox4 = [k for k in ('I2','Cu','Ag','Br2','Cl2') if E[k] > E['Fe3/Fe2']]
print(ox4); assert ox4 == ['Ag','Br2','Cl2']
# 誤った主張： 2*0.80 - 0.34 = 1.26
assert round(2*0.80-0.34,2) == 1.26
print('THE-ELC-047 OK')
