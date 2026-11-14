# THE-RDX-038 検証：FeS の純度と H2S の妨害
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
from sympy import Rational as Q, N
check_item_equations('/home/claude/chemx/todai/out/items/THE-RDX-038/item.json')
M = Q('55.8') + Q('32.1')
mn = Q('0.0200')*Q('50.0')/1000; fe = 5*mn; pur = fe*M/Q('0.500')*100
print('(1)', N(pur)); assert round(float(pur), 1) == 87.9
extra = Q('0.0200')*Q('20.0')/1000; h2s = extra*Q(5, 2)
print('(2)', h2s, h2s/fe*100); assert h2s == Q(1, 1000) and h2s/fe*100 == 20
app = 5*Q('0.0200')*Q('70.0')/1000*M/Q('0.500')*100
print('apparent', N(app)); assert round(float(app)) == 123
print('OK')
