# 検証スクリプト THE-ELC-051（グループ g07）  python3 -I この.py
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

it = check_item('/home/claude/chem6/outh/items/THE-ELC-051/item.json', 4, 5)

import math
E0 = {'Fe': -0.88, 'Ni': 0.49, 'H2': -0.83, 'O2': 0.40}
mn = {'Fe': (2,2), 'Ni': (1,1), 'H2': (2,2), 'O2': (4,4)}
def E(k, oh): m, n = mn[k][1], mn[k][0]; return E0[k] - 0.059*m/n*math.log10(oh)
for oh in [1e-14, 1e-7, 1e-2, 1, 5]:
    assert abs((E('Ni',oh)-E('Fe',oh)) - 1.37) < 1e-9
    assert abs((E('H2',oh)-E('Fe',oh)) - 0.05) < 1e-9
    assert abs((E('Ni',oh)-E('O2',oh)) - 0.09) < 1e-9
print('pH0 換算の水素電位', E('H2', 1e-14))
assert abs(E('H2',1e-14)) < 0.01
F = 96500; ne = 5.00*2.00*3600/F
fe = ne - 2*0.0300; ni = ne - 4*0.0100
print(ne, fe*F/3600, ni*F/3600, fe/ne)
assert round(fe*F/3600,2) == 8.39 and round(ni*F/3600,2) == 8.93 and round(100*fe/ne,1) == 83.9

print('ALL OK')
