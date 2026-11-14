# 検証スクリプト THE-RDX-034（グループ g12）  python3 -I この.py
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


it = check_item('/home/claude/chem6/outh/items/THE-RDX-034/item.json',
                extra_eqs=['Ag + I- -> AgI + e-', '2H+ + 2e- -> H2'])
import math
lg = {'Cl': math.log10(1.8)-10, 'Br': math.log10(5.2)-13, 'I': math.log10(8.3)-17}
lg_given = {'Cl': 0.26-10, 'Br': 0.72-13, 'I': 0.92-17}
for X in lg:
    E = 0.80 + 0.059*lg[X]; Eg = 0.80 + 0.059*lg_given[X]
    print(X, round(E, 4), round(Eg, 4))
    assert round(E, 2) == round(Eg, 2)
E = {X: 0.80 + 0.059*lg_given[X] for X in lg}
assert (round(E['Cl'], 2), round(E['Br'], 2), round(E['I'], 2)) == (0.23, 0.08, -0.15)
dissolve = [X for X in E if E[X] < 0.0]; print('溶ける:', dissolve); assert dissolve == ['I']
# (3) HBr の下限
E0Br = 0.80 + 0.059*lg['Br']; E0Brg = E['Br']
c1 = 10**(E0Br/0.118); c2 = 10**(E0Brg/0.118); print('c下限', c1, c2); assert round(c1) == 4 and round(c2) == 4
assert abs(E0Brg - 0.0755) < 1e-4 and round(E0Brg/0.118, 2) == 0.64
# (4) KI pH 7
EH = 0.059*(-7.0); print('EH', EH); assert round(EH, 2) == -0.41 and E['I'] > EH
print('ALL OK')
