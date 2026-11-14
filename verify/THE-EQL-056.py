# 検証スクリプト THE-EQL-056  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-EQL-056/item.json'
import json, re
from sympy import symbols, solve, simplify
it = json.load(open(ITEM, encoding='utf-8'))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 5, (n_items, n_ans)
eqs = [b for k in ('problem', 'answer', 'solution') for b in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[k])
       if len(re.findall(r'->|<=>|<-', b)) >= 1 and not re.search(r'boxed|\\text|\\mathrm|\$', b)]
print('検査対象の式', eqs); assert eqs == []
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
E, ES, EI, ESI, E0, S, I, K, Ki, k3 = symbols('E ES EI ESI E0 S I K Ki k3', positive=True)
V = k3*E0
# 型①
sol = solve([E*S - K*ES, E*I - Ki*EI, E + ES + EI - E0], [E, ES, EI], dict=True)[0]
assert simplify(k3*sol[ES] - V*S/(K*(1 + I/Ki) + S)) == 0
# 型②
sol = solve([E*S - K*ES, E*I - Ki*EI, ES*I - Ki*ESI, E + ES + EI + ESI - E0], [E, ES, EI, ESI], dict=True)[0]
assert simplify(k3*sol[ES] - V/(1 + I/Ki)*S/(K + S)) == 0
# データの生成
Vv, Kv, Iv = 6.0e-6, 2.0e-4, 2.0e-5
f0 = lambda s: Vv*s/(Kv + s)
fX = lambda s: Vv*s/(Kv*(1 + Iv/1.0e-5) + s)
fY = lambda s: Vv/(1 + Iv/4.0e-5)*s/(Kv + s)
data = {'none': (2.00e-6, 5.00e-6), 'X': (8.57e-7, 3.75e-6), 'Y': (1.33e-6, 3.33e-6)}
for key, f in (('none', f0), ('X', fX), ('Y', fY)):
    for s, given in zip((1e-4, 1e-3), data[key]):
        assert sig(f(s), 3) == given, (key, s, f(s))
def fit(v1, v2, s1=1e-4, s2=1e-3):
    slope = (1/v1 - 1/v2)/(1/s1 - 1/s2); icpt = 1/v2 - slope/s2
    Vp = 1/icpt; return Vp, slope*Vp
res = {k: fit(*v) for k, v in data.items()}; print(res)
assert [sig(x, 2) for x in res['none']] == [6.0e-6, 2.0e-4]
assert [sig(x, 2) for x in res['X']] == [6.0e-6, 6.0e-4]
assert [sig(x, 2) for x in res['Y']] == [4.0e-6, 2.0e-4]
aX = res['X'][1]/res['none'][1]; aY = res['none'][0]/res['Y'][0]; print('alpha', aX, aY)
assert sig(Iv/(aX - 1), 2) == 1.0e-5 and sig(Iv/(aY - 1), 2) == 4.0e-5
s = 1e-2; print('X', fX(s)/f0(s), 'Y', fY(s)/f0(s))
assert round(fX(s)/f0(s)*100) == 96 and round(fY(s)/f0(s)*100) == 67
print('全検証OK:', ITEM)
