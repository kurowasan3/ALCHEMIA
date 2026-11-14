# 検証スクリプト THE-EQL-089（グループ g12）  python3 -I この.py
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


it = check_item('/home/claude/chem6/outh/items/THE-EQL-089/item.json')
from sympy import Rational as Q, symbols, solve
p = Q(2, 3)
y = symbols('y', positive=True)
nh = Q(1,2) + y; nx = 2 - 2*y
HD = nh*2*p*(1-p)
assert 2*p*(1-p) == Q(4,9)
sol = solve(HD/y - Q(4,3), y); print('y =', sol); assert sol == [Q(1,4)]
yv = sol[0]; nhv = nh.subs(y, yv); nxv = nx.subs(y, yv)
# 原子の保存
assert nxv + 2*yv == 2 and nxv + 2*nhv == 3
H2 = nhv*p**2; HDv = nhv*2*p*(1-p); D2 = nhv*(1-p)**2; HI = nxv*p; DI = nxv*(1-p); I2 = yv
vals = dict(H2=H2, HD=HDv, D2=D2, HI=HI, DI=DI, I2=I2)
print({k: float(v) for k, v in vals.items()})
assert (H2, HDv, D2, HI, DI, I2) == (Q(1,3), Q(1,3), Q(1,12), 1, Q(1,2), Q(1,4))
# 原子の総数（H, D, I）
assert 2*H2 + HDv + HI == 2 and 2*D2 + HDv + DI == 1 and HI + DI + 2*I2 == 2
K = HI**2/(H2*I2); print('K =', K); assert K == 12
assert DI**2/(D2*I2) == 12 and nxv**2/(nhv*yv) == 12
assert HDv**2/(H2*D2) == 4
Kp = HI*DI/(HDv*I2); print("K' =", Kp); assert Kp == 6
# 解の一意性：HD/I2 は y について単調減少
import sympy
assert sympy.diff(HD/y, y).subs(y, Q(1,4)) < 0
# 有効数字 2 桁の表記
assert ['%.2g' % float(v) for v in (H2, HDv, D2, HI, DI, I2)] == ['0.33', '0.33', '0.083', '1', '0.5', '0.25']
print('ALL OK')
