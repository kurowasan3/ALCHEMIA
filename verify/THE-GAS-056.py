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

it = check_item('/home/claude/chemx/kyodai/outh/items/THE-GAS-056/item.json')
RT=8.3e3*77
def run(PHe_after=1.25e4):
    VI = 0.100*(2.00e4/PHe_after-1)
    N1 = (2.88e3*0.100 - 1.00e3*(0.100+VI))/RT
    N2 = (7.08e3*0.100 + 1.00e3*VI - 4.00e3*(0.100+VI))/RT + N1
    N0 = (4.00e3-1.00e3)/(4.00e3/N2 - 1.00e3/N1)
    K = 1/((1.00e3/N1 - 1.00e3/N0)*N0)
    S = N0*6.0e23*1.6e-19/0.050
    return VI,N1,N2,N0,K,S
VI,N1,N2,N0,K,S = run()
print(VI,N1,N2,N0,K,S)
assert abs(VI-0.060)<1e-9 and sig(N1)==2.0e-4 and sig(N2)==4.0e-4 and sig(N0)==6.0e-4 and sig(K)==5.0e-4 and sig(S)==1.2e3
for P,N in [(1.00e3,N1),(4.00e3,N2)]:
    assert abs(N0*K*P/(1+K*P)-N)/N < 0.01
# ヘリウムがわずかに吸着 → 圧力がさらに下がる
VIb,_,_,N0b,_,Sb = run(1.24e4)
assert VIb>VI and Sb<S
print('He吸着時', VIb, Sb)
import xml.dom.minidom; xml.dom.minidom.parse('/home/claude/chemx/kyodai/outh/items/THE-GAS-056/problem.svg')
print('OK')
