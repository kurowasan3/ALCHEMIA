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

it = check_item('/home/claude/chemx/kyodai/outh/items/THE-ION-097/item.json')
from scipy.optimize import brentq
Ka,S,Cc,K1,K2,Kw=1.0e-5,2.0e-6,1.0e-5,4.5e-7,4.7e-11,1e-14
K=Ka*S
def solve(C):
    f=lambda lh: (lambda h: C+h-(K/h+K1*Cc/h+2*K1*K2*Cc/h**2+Kw/h))(10**lh)
    h=10**brentq(f,-14,-1); ppt=C-K/h-S
    if ppt>=0: return True,-math.log10(h),ppt/C
    g=lambda lh: (lambda h: C+h-(C*Ka/(Ka+h)+K1*Cc/h+2*K1*K2*Cc/h**2+Kw/h))(10**lh)
    h=10**brentq(g,-14,-1); return False,-math.log10(h),0
for C,ph,x in [(1e-3,7.6,18),(1e-4,6.6,16)]:
    pp,p,xx=solve(C); h=(K+K1*Cc)/C; xa=K1*Cc/(K+K1*Cc)-S/C
    print(C,p,xx,-math.log10(h),xa)
    assert pp and round(p,1)==ph and round(-math.log10(h),1)==ph and round(xa*100)==x
    assert abs(xx-xa)<0.005
pp,p,_=solve(1e-5); ha=math.sqrt(K1*Cc/(1+1e-5/Ka))
assert (not pp) and round(p,1)==5.8 and round(-math.log10(ha),1)==5.8
assert K1*Cc/(K+K1*Cc)-S/1e-5<0
assert round(-math.log10(math.sqrt(K1*Cc)),1)==5.7
assert abs(K2/2.45e-8-1.9e-3)<1e-4
assert abs(math.log10(2.45)-0.39)<0.005 and abs(math.log10(1.5)-0.18)<0.005
print('OK')
