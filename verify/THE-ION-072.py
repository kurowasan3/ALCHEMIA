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

it = check_item('/home/claude/chem6/out/items/THE-ION-072/item.json')
from sympy import symbols, solve, sqrt, Rational as Q, nsolve
K1, K2 = Q(1,10**7), Q(1,10**14)
Ksp = {'FeS': Q(5,10**18), 'ZnS': Q(2,10**24), 'CuS': Q(6,10**36)}
K = {k: v/(K1*K2) for k, v in Ksp.items()}
print({k: float(v) for k, v in K.items()})
assert [sig(float(K[k])) for k in ('FeS','ZnS','CuS')] == [5.0e3, 2.0e-3, 6.0e-15]
def full(cHCl, n=Q(1,100)):
    Qr = n*n/(cHCl-2*n)**2
    return {k: Qr < K[k] for k in K}, float(Qr)
r1, q1 = full(Q(1)); r2, q2 = full(Q(1,10)); print(r1, q1, r2, q2)
assert r1 == {'FeS':True,'ZnS':True,'CuS':False} and r2 == {'FeS':True,'ZnS':False,'CuS':False}
s = symbols('s', positive=True)
sz = [x for x in solve(s**2 - K['ZnS']*(Q(1,10)-2*s)**2, s) if x < Q(1,20)][0]
print('ZnS dissolved', float(sz)); assert sig(float(sz)) == 4.1e-3 and float(sz) < 0.01
sc = nsolve(s**2 - float(K['CuS'])*(1-2*s)**2, s, 1e-7); print('CuS', sc); assert sig(float(sc)) == 7.7e-8
# H2S が飽和（0.10 mol/L）しないこと
assert 0.01 < 0.10
for e in ['3CuS + 8H+ + 2NO3- -> 3Cu^2+ + 3S + 2NO + 4H2O', 'NO3- + 4H+ + 3e- -> NO + 2H2O', 'S^2- -> S + 2e-']:
    assert check_equation(e) is None
print('THE-ION-072 OK')
