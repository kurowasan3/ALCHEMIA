# THE-RDX-042 検証：FeS + Fe2O3 試料の組成
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
it = check_item_equations('/home/claude/chemx/todai/outh/items/THE-RDX-042/item.json')
for eq in ['Fe2O3 + 3H2SO4 -> Fe2(SO4)3 + 3H2O', 'FeS + H2SO4 -> FeSO4 + H2S', 'CuSO4 + H2S -> CuS + H2SO4']:
    assert check_equation(eq) is None, eq
Fe, S, O, Cu = Q('55.8'), Q('32.1'), Q('16.0'), Q('63.5')
x, y = symbols('x y')
sol = solve([x - y - Q('0.478')/(Cu + S), x + 2*y - 5*Q('0.0200')*Q('0.0800')], [x, y])
print(sol); assert sol[x] == Q(6, 1000) and sol[y] == Q(1, 1000)
pFeS = sol[x]*(Fe + S)*100; pFe2O3 = sol[y]*(2*Fe + 3*O)*100
print(N(pFeS), N(pFe2O3)); assert round(float(pFeS), 1) == 52.7 and round(float(pFe2O3), 1) == 16.0
i = Q(8, 1000)*(Fe + S)*100; ii = Q(5, 1000)*(Fe + S)*100
print(N(i), N(ii)); assert round(float(i), 1) == 70.3 and round(float(ii), 1) == 44.0
x2 = 5*Q('0.0200')*Q('0.0300')/3; print('x2', x2, N(x2*(Fe + S)), N(x2*(2*Fe + 3*O)))
assert x2 == Q(1, 1000) and round(float(x2*(Fe + S)), 4) == 0.0879
print('OK')
