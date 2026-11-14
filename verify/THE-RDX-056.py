# 検証スクリプト（グループ g11）  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-RDX-056/item.json'
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
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
eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
assert check_equation('2XeF2 + 2H2O -> 2Xe + 3HF + O2') is not None  # 検査器が誤りを検出できること
def sig(x, n):
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
from scipy.optimize import brentq
E1, E2, k = 0.68, 1.44, 0.059
logK = (E2-E1)/k; K = 10**logK; print(logK, K); assert sig(K, 2) == 7.6e12
assert round(10**0.88, 1) == 7.6 and abs(10**0.44 - 2.75) < 0.01 and abs(10**0.356 - 2.27) < 0.005
def E(v):  # 物質収支と平衡定数を連立して数値的に解く
    Fe, Ce = 2.0, 0.1*v
    x = brentq(lambda x: x*x - K*(Fe-x)*(Ce-x), 0, min(Fe, Ce)*(1-1e-14))
    e_fe = E1 + k*math.log10(x/(Fe-x))
    return e_fe
for v, ans in ((5.0, 0.65), (10.0, 0.68), (18.0, 0.74), (19.98, 0.86), (20.0, 1.06)):
    e = E(v); print(v, e); assert round(e, 2) == ans
def E_after(v):
    Fe, Ce = 2.0, 0.1*v
    x = brentq(lambda x: x*x - K*(Fe-x)*(Ce-x), 0, Fe*(1-1e-14))
    return E2 + k*math.log10((Ce-x)/x)
for v, ans in ((20.02, 1.26), (22.0, 1.38), (40.0, 1.44)):
    e = E_after(v); print(v, e); assert round(e, 2) == ans
fe2 = 0.0500/math.sqrt(K); print('Fe2+ eq', fe2); assert sig(fe2, 2) == 1.8e-8
r = 10**((0.76-E1)/k); frac = r/(1+r)*100; print(frac); assert sig(frac, 3) == 95.8
# MnO4- の当量点：平衡組成を数値的に解いて式と比べる
for H, ans in ((1.0, 1.37), (0.10, 1.29)):
    formula = (E1 + 5*1.51)/6 + 8*k/6*math.log10(H)
    # 当量点：Fe 2.00 mmol, MnO4- 0.400 mmol, 体積は相殺されるので比で考える（濃度は 40 mL で割る）
    V = 0.040; Fe0, Mn0 = 2.0e-3/V*1e-0/1, 0.4e-3/V
    Fe0 = 2.0e-3/ (0.020+0.020); Mn0 = Fe0/5
    logKMn = 5*(1.51-E1)/k
    # x = Mn2+ 生成濃度：K' = [Mn2+][Fe3+]^5/([MnO4-][Fe2+]^5[H+]^8)
    # 未反応の MnO4- 濃度 z を対数で解く（z は非常に小さいので log z を未知数にする）
    f = lambda lz: math.log10(Mn0-10**lz) + 5*math.log10(5*(Mn0-10**lz)) - lz - 5*math.log10(5*10**lz) - 8*math.log10(H) - logKMn
    lz = brentq(f, math.log10(Mn0)-40, math.log10(Mn0)-1e-9, xtol=1e-14)
    z = 10**lz
    e_num = E1 + k*math.log10((Mn0-z)/z)
    print(H, formula, e_num); assert abs(formula - e_num) < 1e-6 and round(formula, 2) == ans
assert check_equation('MnO4- + 5Fe^2+ + 8H+ -> Mn^2+ + 5Fe^3+ + 4H2O') is None
print('OK')
