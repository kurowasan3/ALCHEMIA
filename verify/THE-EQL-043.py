# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-EQL-043/item.json'
from sympy import *
import sys, json, re
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
print('反応式', len(eqs), '本を確認')

def sig(x, n=2):
    x = float(x)
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def chk(name, got, want, n=2):
    assert sig(got, n) == sig(want, n), (name, float(got), want)
    print('OK', name, sig(got, n))
K1, K2 = 4.0e3, 1.0e5
c = 8.3/(8.3e3*1000); assert abs(c - 1e-6) < 1e-15
def solve_sys(nCa, nC):
    """平衡状態を場合分けで求める（一意な整合解を返す）"""
    z = symbols('z', positive=True)
    out = []
    # 両方残る
    pCO2, pCO = K1, (K2*K1)**0.5
    n1, n2 = (pCO2 + pCO/2)*c, pCO/2*c
    if n1 < nCa and n2 < nC: out.append(('both', pCO2, pCO, nCa-n1, nC-n2))
    # CaCO3 消失，黒鉛あり
    zz = [s for s in solve(Eq(z**2/K2 + z/2, nCa/c), z)][0]
    pCO2, pCO = float(zz**2/K2), float(zz)
    if pCO2 < K1 and pCO/2*c < nC: out.append(('noCa', pCO2, pCO, 0, nC - pCO/2*c))
    # 黒鉛消失，CaCO3 あり
    pCO2, pCO = K1, 2*nC/c
    if pCO**2/pCO2 < K2 and (pCO2 + pCO/2)*c < nCa: out.append(('noC', pCO2, pCO, nCa - (pCO2+pCO/2)*c, 0))
    assert len(out) == 1, out
    return out[0]
r = solve_sys(0.050, 0.050); print(r)
chk('(1) CO2', r[1], 4.0e3); chk('(1) CO', r[2], 2.0e4); chk('(1) P', r[1]+r[2], 2.4e4); chk('(1) CaCO3', r[3], 3.6e-2); chk('(1) C', r[4], 4.0e-2)
chk('(2) CaCO3 min', (K1 + (K1*K2)**0.5/2)*c, 1.4e-2); chk('(2) C min', (K1*K2)**0.5/2*c, 1.0e-2)
r = solve_sys(6.0e-3, 0.050); print(r); assert r[0] == 'noCa'
chk('(3) CO2', r[1], 1.0e3); chk('(3) CO', r[2], 1.0e4); chk('(3) P', r[1]+r[2], 1.1e4)
r = solve_sys(0.050, 5.0e-3); print(r); assert r[0] == 'noC'
chk('(4) CO2', r[1], 4.0e3); chk('(4) CO', r[2], 1.0e4); chk('(4) P', r[1]+r[2], 1.4e4); chk('(4) CaCO3', r[3], 4.1e-2)
print('THE-EQL-043 すべてOK')
