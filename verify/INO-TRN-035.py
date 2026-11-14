# INO-TRN-035 検証：CuO/Cu2O 混合物のヨウ素滴定
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
check_item_equations('/home/claude/chemx/todai/out/items/INO-TRN-035/item.json')
assert check_equation('2Cu + O2 -> 2CuO') is None
Cu, O = Q('63.5'), Q('16.0')
s = Q('0.100')*Q('0.0300'); i2 = s/2; cuo = 2*i2
m_cuo = cuo*(Cu + O); cu2o = (Q('0.310') - m_cuo)/(2*Cu + O)
cu = cuo + 2*cu2o; pct = cu*Cu/Q('0.310')*100
print(cuo, N(m_cuo), N(cu2o), N(pct))
assert cuo == Q(3, 1000) and cu2o == Q(5, 10000) and round(float(pct), 1) == 81.9
print('frac', 2*cu2o/cu); assert 2*cu2o/cu == Q(1, 4)
print('OK')
