# 検証スクリプト THE-EQL-081（グループ g12）  python3 -I この.py
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq): return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq): return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1: return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

FIELDS = set("""原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸""".split())

def check_item(path, extra_eqs=()):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            e = clean_ce(b)
            if not should_check(b): continue
            r = check_equation(e)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    for e in extra_eqs:
        r = check_equation(e); assert r is None, (e, r); n_eq += 1; print('式OK(追加):', e)
    print('検査した式の数:', n_eq)
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    assert it['field'] in FIELDS, it['field']
    for f in it['format'].split(','): assert f in ('qa','equation','compare','condition','pathway','calc','structure'), f
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it


it = check_item('/home/claude/chem6/out/items/THE-EQL-081/item.json')
from sympy import Rational as Q, symbols, solve, simplify
p = symbols('p', positive=True)
K = simplify((2*p*(1-p))**2/(p**2*(1-p)**2)); print('K =', K); assert K == 4
# 仕込み H2 1.0, D2 4.0
x = symbols('x')
sols = [s for s in solve(4*x**2 - 4*(1-x)*(4-x), x) if 0 < s < 1]
print('x =', sols); assert sols == [Q(4,5)]
xv = sols[0]; H2, HD, D2 = 1-xv, 2*xv, 4-xv
pv = Q(2,10); N = 5
assert (H2, HD, D2) == (N*pv**2, N*2*pv*(1-pv), N*(1-pv)**2)
print('H2,HD,D2 =', float(H2), float(HD), float(D2)); assert (float(H2), float(HD), float(D2)) == (0.2, 1.6, 3.2)
M = (1*2 + 4*4)/Q(5); print('平均分子量', float(M)); assert M == Q(18,5)
Mafter = (H2*2 + HD*3 + D2*4)/(H2+HD+D2); assert Mafter == M
r = [H2/H2, HD/H2, D2/H2]; print('ピーク比', r); assert r == [1, 8, 16]
print('ALL OK')
