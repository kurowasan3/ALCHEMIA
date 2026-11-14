# THE-EQL-103 検証：定常状態近似による HI 生成速度
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
from sympy import symbols, solve, Rational as Q, simplify, exp, diff, log, N
check_item_equations('/home/claude/chemx/todai/outh/items/THE-EQL-103/item.json')
kd, kr, k3, H, I2, I_2 = symbols('k_d k_r k_3 H I2 Isq', positive=True)
Isq = solve(2*kd*I2 - 2*kr*I_2 - 2*k3*H*I_2, I_2)[0]
v = simplify(2*k3*H*Isq); print('v =', v)
assert simplify(v - 2*kd*k3*H*I2/(kr + k3*H)) == 0
vm, K = symbols('vm K', positive=True)
sol = solve([vm*Q('0.050')/(K + Q('0.050')) - Q('2.0e-6'), vm*Q('0.20')/(K + Q('0.20')) - Q('5.0e-6')], [vm, K], dict=True)[0]
print(sol); assert sol[K] == Q(1, 5) and sol[vm] == Q(1, 10**5)
print('kd', sol[vm]/(2*Q('1.0e-3'))); assert sol[vm]/(2*Q('1.0e-3')) == Q(5, 1000)
v1 = sol[vm]*1/(sol[K] + 1); print('v(1.0)', N(v1)); assert round(float(v1)*1e6, 1) == 8.3
# 見かけの活性化エネルギー（数値微分）
R = 8.314e-3; T = symbols('T', positive=True)
Ed, Er, E3 = 155, 5, 20
def k(E, A=1.0): return A*exp(-E/(R*T))
for Hval, A3, label in [(1e-12, 1.0, 'low'), (1e12, 1.0, 'high')]:
    vv = 2*k(Ed)*k(E3, A3)*Hval/(k(Er) + k(E3, A3)*Hval)
    Ea = -R*diff(log(vv), T)*(-T**2)
    val = float(Ea.subs(T, 600)); print(label, round(val, 2))
    assert abs(val - (170 if label == 'low' else 155)) < 0.5
ratio = k(E3)/k(Er); assert float(ratio.subs(T, 700)) > float(ratio.subs(T, 600))
print('OK')
