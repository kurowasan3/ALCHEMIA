# INO-HAL-017 検証：水素化物の結晶中の水素密度と加水分解
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
check_item_equations('/home/claude/chemx/todai/out/items/INO-HAL-017/item.json')
NA = Q('6.0e23'); mH = 4*Q('1.0')/NA
def vol(d_nm): a = 2*Q(d_nm)*Q('1e-7'); return a**3/mH
vL, vN, vliq = vol('0.204'), vol('0.244'), 1/Q('0.071')
print('LiH', N(vL), 'NaH', N(vN), 'liq', N(vliq))
assert round(float(vL)) == 10 and round(float(vN)) == 17 and round(float(vliq)) == 14
assert vL < vliq < vN
n = 1/Q('7.9'); mH2 = n*2; print('H2 g', N(mH2), 'ratio', mH2/(n*1))
assert round(float(mH2), 2) == 0.25 and mH2/(n*1) == 2
print('OK')
