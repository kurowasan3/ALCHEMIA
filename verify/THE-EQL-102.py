# THE-EQL-102 検証：一定圧力シリンダーと膜でつながった N2O4/NO2 平衡
# --- 共通：item.json の \ce{} 反応式のつり合い検査（単独で動く） ---
import json, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
def _ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k, d = j + 4, 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k
def check_item_equations(path):
    it = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for ce in _ce_blocks(it[key]):
            arrows = re.findall(r'<=>|->|<-', ce)
            if len(arrows) != 1 or '\\boxed' in ce or '\\text' in ce or '\\mathrm' in ce:
                continue
            r = check_equation(ce)
            assert r is None, (key, ce, r)
            n += 1
    print('反応式', n, '本 OK')
    return it
from sympy import Rational as Q, symbols, solve, N
check_item_equations('/home/claude/chemx/todai/outh/items/THE-EQL-102/item.json')
x = symbols('x', positive=True)
RT = Q(25, 10)*10**6; V = 10; Kp = 10**5
def closed(n0, Kp=Kp, RT=RT):
    P0 = n0*RT/V
    s = [r for r in solve(4*x**2 - Kp*(P0 - x), x) if 0 < r < P0][0]
    return 2*s, P0 - s
pn, pd = closed(Q('0.40')); print('(1) closed p_NO2', N(pn)); assert pn > Q(5, 10)*10**5
p0 = Q(5, 10)*10**5; pdf = p0**2/Kp
nbox_N = (2*pdf + p0)*V/RT; ncyl = 2*Q('0.40') - nbox_N; Vcyl = ncyl*RT/p0
print('(1)', pdf, ncyl, Vcyl); assert pdf == 25000 and ncyl == Q(2, 5) and Vcyl == 20
pn2, _ = closed(Q('0.10')); print('(2)', N(pn2)); assert pn2 < p0 and round(float(pn2)/1e4, 1) == 3.1
n_th = (pdf + p0/2)*V/RT; print('(3)', n_th); assert n_th == Q(1, 5)
print('(4)', N(pn)); assert round(float(pn)/1e4, 1) == 7.8
# (5) 温度上昇：Kp と RT がともに大きくなる例
Kp2, RT2 = Q(12, 10)*10**5, RT*Q(102, 100)
nbox2 = (2*p0**2/Kp2 + p0)*V/RT2; print('(5) box N', N(nbox_N), '->', N(nbox2)); assert nbox2 < nbox_N
print('OK')
