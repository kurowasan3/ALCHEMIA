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

it = check_item('/home/claude/chem6/outh/items/THE-THC-041/item.json')
from sympy import Rational, nsimplify
for e in ['2C + O2 -> 2CO','C + O2 -> CO2','2Ni + O2 -> 2NiO','2Fe + O2 -> 2FeO','4/3Cr + O2 -> 2/3Cr2O3','4/3Al + O2 -> 2/3Al2O3','C + CO2 -> 2CO','2NiO + C -> 2Ni + CO2','2/3Cr2O3 + 4/3Al -> 4/3Cr + 2/3Al2O3']:
    assert check_equation(e) is None, e
# 熱力学データから表の値を再計算
fH = dict(C=0,O2=0,CO=-110.5,CO2=-393.5,Ni=0,NiO=-239.7,Fe=0,FeO=-272.0,Cr=0,Cr2O3=-1139.7,Al=0,Al2O3=-1675.7)
S0 = dict(C=5.7,O2=205.2,CO=197.7,CO2=213.8,Ni=29.9,NiO=38.0,Fe=27.3,FeO=60.75,Cr=23.8,Cr2O3=81.2,Al=28.3,Al2O3=50.9)
rx = {1:({'C':2,'O2':1},{'CO':2}),2:({'C':1,'O2':1},{'CO2':1}),3:({'Ni':2,'O2':1},{'NiO':2}),4:({'Fe':2,'O2':1},{'FeO':2}),5:({'Cr':Rational(4,3),'O2':1},{'Cr2O3':Rational(2,3)}),6:({'Al':Rational(4,3),'O2':1},{'Al2O3':Rational(2,3)})}
table = {1:(-221,179),2:(-394,3),3:(-479,-189),4:(-544,-138),5:(-760,-183),6:(-1117,-209)}
for k,(L,R) in rx.items():
    h = sum(fH[s]*c for s,c in R.items())-sum(fH[s]*c for s,c in L.items())
    sS = sum(S0[s]*c for s,c in R.items())-sum(S0[s]*c for s,c in L.items())
    print(k, round(float(h),1), round(float(sS),1), table[k])
    assert abs(float(h)-table[k][0]) <= 1.0 and abs(float(sS)-table[k][1]) <= 1.0
T = lambda a, b: (table[a][0]-table[b][0]) / ((table[a][1]-table[b][1])/1000)
t12 = T(1,2); print('C+CO2->2CO', t12); assert round(t12) == 983
res = {}
for m in (3,4,5):
    tco, tco2 = T(1,m), T(2,m)
    best = min((tco,'CO'),(tco2,'CO2')); res[m] = best
    print(m, round(tco,1), round(tco2,1), best)
assert f'{res[3][0]:.3g}' == '443' and res[3][1] == 'CO2'
assert f'{res[4][0]:.3g}' == '1.02e+03' and res[4][1] == 'CO'
assert f'{res[5][0]:.3g}' == '1.49e+03' and res[5][1] == 'CO'
# 還元反応の ΔS>0（C 経路）・ΔH>0 で T 下限になっていることの確認
for m in (3,4,5):
    for c in (1,2):
        assert table[c][0]-table[m][0] > 0 and table[c][1]-table[m][1] > 0
dH = table[6][0]-table[5][0]; dS = table[6][1]-table[5][1]
g300 = dH - 300*dS/1000; g2000 = dH - 2000*dS/1000; print('thermite', dH, dS, g300, g2000)
assert round(g300) == -349 and round(g2000) == -305
print('THE-THC-041 OK')
