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

it = check_item('/home/claude/chemx/kyodai/outh/items/THE-ION-096/item.json')
from scipy.optimize import brentq
K1,K2,K3,Kw=7.1e-3,6.3e-8,4.5e-13,1.0e-14
def exact(C, na):
    def f(lh):
        h=10**lh; D=h**3+K1*h*h+K1*K2*h+K1*K2*K3
        a1=K1*h*h/D; a2=K1*K2*h/D; a3=K1*K2*K3/D
        return na*C + h - (C*(a1+2*a2+3*a3) + Kw/h)
    return -brentq(f,-14,0)
ph = lambda h: -math.log10(h)
p1 = ph(math.sqrt(K2*(K3*0.10+Kw)/(K2+0.10))); p2 = ph(math.sqrt(K2*(K3*1e-3+Kw)/(K2+1e-3)))
p4 = ph(math.sqrt(K1*(K2*1e-3+Kw)/(K1+1e-3)))
print(p1, exact(0.10,2), p2, exact(1e-3,2), p4, exact(1e-3,1))
assert round(p1,1)==9.7 and round(p2,1)==9.1 and round(p4,1)==5.1
assert abs(p1-exact(0.10,2))<0.01 and abs(p2-exact(1e-3,2))<0.01 and abs(p4-exact(1e-3,1))<0.01
assert round((ph(K2)+ph(K3))/2,1)==9.8 and round((ph(K1)+ph(K2))/2,1)==4.7
for v,l in [(3.47,0.54),(6.58,0.82),(5.52,0.74),(6.3,0.80),(4.5,0.65),(7.1,0.85)]:
    assert abs(math.log10(v)-l)<0.006, v
print('Kw/K3 =', Kw/K3)
print('OK')
