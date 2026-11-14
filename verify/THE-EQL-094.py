# THE-EQL-094 検証：膜でつないだ N2O4/NO2 平衡箱
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
from sympy import Rational as Q, sqrt, N
check_item_equations('/home/claude/chemx/todai/out/items/THE-EQL-094/item.json')
Kc = Q('0.40')
c1 = Q('0.50')/5; c2 = sqrt(Kc*c1)
print('[N2O4]', c1, '[NO2]', N(c2)); assert c1 == Q(1,10) and c2 == Q(1,5)
# 原子収支：加えた N2O4 0.50 mol と取り出した NO2 1.0 mol
assert 2*Q('0.50') == 1*Q('1.0') and 4*Q('0.50') == 2*Q('1.0')
VQ = Q('1.0')/c2; print('VQ', VQ); assert VQ == 5
c1b = Q('0.50')/20; c2b = sqrt(Kc*c1b); VQb = 1/c2b
print('(5)', N(c2b), N(VQb)); assert VQb == 10
print('OK')
