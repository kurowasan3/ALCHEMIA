# 検証スクリプト INO-PPT-034  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/INO-PPT-034/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    return re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', s)
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq
def should_check(eq):
    # 矢印が1つの完全な式はすべて（validate ツールが一般式とみなして飛ばす Zn・Sn なども含めて）検査する．
    # 一般式・同位体の式は矢印に空ラベル [\\text{}] を付けてあるので除く．
    if re.search(r'boxed|□|\?|？|…|\\square', eq): return False
    if re.search(r'(->|<=>|<-)\[\\text\{\}\]', eq): return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1: return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None
def check_all_equations():
    eqs = sorted({b for k in ('problem','answer','solution') for b in ce_blocks(it[k]) if should_check(b)})
    for e in eqs:
        r = check_equation(clean_ce(e)); print('式', e, '→', r or 'OK'); assert r is None, e
    return eqs
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items+1)), nums
    print('小問数', n_items)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
EQS = check_all_equations()

# ---- 個別の検証 ----

from sympy import Rational as R, sqrt, symbols, solve, log, N
gen = [b for k in ('problem','solution') for b in ce_blocks(it[k]) if '[\\text{}]' in b]
for g in set(gen):
    r = check_equation(clean_ce(g)); print('一般式', g, r or 'OK'); assert r is None
# α を EDTA の pKa から計算して表と比べる
pKa = [2.0, 2.67, 6.16, 10.26]
Ks = [10**-p for p in pKa]
def alpha(pH):
    h = 10**-pH; K1,K2,K3,K4 = Ks
    num = K1*K2*K3*K4
    return num/(h**4 + K1*h**3 + K1*K2*h**2 + K1*K2*K3*h + num)
table = {2:3.7e-14, 4:3.6e-9, 6:2.2e-5, 8:5.4e-3, 10:0.35, 12:0.98}
for p, a in table.items():
    print(p, alpha(p), a); assert abs(alpha(p)/a - 1) < 0.1
KCa, KFe = 5.0e10, 1.3e25
C = 5.00e-3
# (1) 同値性：[M]/C = 1/sqrt(K'C)
for Kp in (1e8, 2e8, 1e10): assert abs(((C/Kp)**0.5)/C - 1/(Kp*C)**0.5) < 1e-15
need = 1e6/C; assert abs(need - 2.0e8) < 1; assert abs(need/KCa - 4.0e-3) < 1e-12
ok = [p for p in sorted(table) if table[p]*KCa*C >= 1e6]; print('ok pH', ok); assert ok[0] == 8
f6 = 1/(table[6]*KCa*C)**0.5; print('pH6', f6); assert sig(f6*100,2) == 1.3
f10 = 1/(table[10]*KCa*C)**0.5; print('pH10', f10); assert sig(f10*100,2) == 0.011
# 近似の妥当性：厳密解 x^2 = (C-x)/K' と比較（pH 6）
x = symbols('x', positive=True); Kp6 = table[6]*KCa
ex = [s for s in solve(Kp6*x**2 + x - C, x) if s > 0][0]; ap = (C/Kp6)**0.5
print('厳密', float(ex), '近似', ap); assert abs(float(ex)/ap - 1) < 0.02
KpFe = table[2]*KFe; print('KFe', KpFe); assert sig(KpFe,2) == 4.8e11
assert sig(KpFe*C,2) == 2.4e9
Y = (C/KpFe)**0.5; print('Y', Y); assert sig(Y,3) == 1.02e-7
KpCa2 = table[2]*KCa; print('KCa pH2', KpCa2); assert sig(KpCa2,2) == 1.9e-3
ratio = KpCa2*Y; print('ratio', ratio); assert sig(ratio,2) == 1.9e-10
assert sig(table[10]*KCa,3) == 1.75e10
H = 2*5.0e-5/0.0550; print('H', H); assert sig(H,3) == 1.82e-3
import math; pH = -math.log10(H); print('pH', pH); assert round(pH,1) == 2.7
assert sig(table[4]*KCa,2) == 1.8e2

print('全検証OK:', ITEM)
