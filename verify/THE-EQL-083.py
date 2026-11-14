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

it = check_item('/home/claude/chem6/outh/items/THE-EQL-083/item.json')
import math
from sympy import Rational, log, exp, nsimplify
ratio = Rational(16*24, 12); print('k360/k320', ratio); assert ratio == 32
d1 = Rational(1,320) - Rational(1,360); d2 = Rational(1,300) - Rational(1,320)
print(d1, d2); assert d1 == Rational(1,2880) and d2 == Rational(1,4800)
Ea = 5*0.693*8.31*2880; print('Ea', Ea); assert f'{Ea:.3g}' == '8.29e+04'
Ea_exact = 5*math.log(2)*8.31*2880
r = math.exp(Ea_exact/8.31*float(d2)); print('k320/k300', r); assert abs(r-8) < 1e-9
K = -math.log(0.9); k300 = K/128; k320 = K/16
def day_reach(schedule):     # schedule: list of (k, days) ; 最後は 300 K で到達まで
    s = 0; t = 0
    for kk, dd in schedule:
        if s + kk*dd >= K: return t + (K-s)/kk
        s += kk*dd; t += dd
    return t + (K-s)/k300
a = day_reach([(k300,64),(k320,4)]); b = day_reach([(k320,4)])
print(a, b); assert abs(a-100) < 1e-9 and abs(b-100) < 1e-9
# 二次反応でも順序によらない
for C0 in (0.5, 1.0, 3.0):
    kk300, kk320 = 0.001, 0.008
    def run(seq):
        inv = 1/C0
        for kk, dd in seq: inv += kk*dd
        return 1/inv
    A_ = run([(kk300,64),(kk320,4),(kk300,30)]); B_ = run([(kk320,4),(kk300,94)])
    assert abs(A_-B_) < 1e-12
print('THE-EQL-083 OK')
