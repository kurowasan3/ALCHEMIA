import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_all_equations(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts(it):
    np_ = it['problem'].count('\\item')
    na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), (np_, na)
    print('小問数OK:', np_)
import sympy as sp, numpy as np
it = json.load(open('/home/claude/chem3/outh/items/THE-GAS-038/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
P, K, b, R, T = sp.symbols('P K b R T', positive=True)
pM = (-1 + sp.sqrt(1 + 4*K*P))/(2*K)
assert sp.simplify(pM + K*pM**2 - P) == 0
Z = P/(2*P - pM) + b*P/(R*T)
ser = sp.series(Z, P, 0, 2).removeO()
print('low P:', sp.simplify(ser))
assert sp.simplify(ser - (1 - (K - b/(R*T))*P)) == 0
lim = sp.limit(P/(2*P - pM), P, sp.oo); print('limit', lim); assert lim == sp.Rational(1, 2)
Rv, Tv = 8.3e3, 300
bv = 1.6e-8*Rv*Tv; m = bv/6.0e23; Kv = 2.4e-8 + 1.6e-8
print('b', bv, 'm', m, 'K', Kv); assert sig(m, 2) == 6.6e-26 and sig(Kv, 2) == 4.0e-8
pMv = (-1 + np.sqrt(1 + 4*Kv*5.0e7))/(2*Kv); Zv = 5.0e7/(1.0e8 - pMv) + bv*5.0e7/(Rv*Tv)
frac = 2*(5.0e7 - pMv)/(1.0e8 - pMv)
print('pM', pMv, 'Z', Zv, 'frac', frac); assert abs(pMv - 2.5e7) < 1 and sig(Zv, 2) == 1.5 and abs(frac - 2/3) < 1e-9
def Zf(Pa, Kk, TT):
    p = (-1 + np.sqrt(1 + 4*Kk*Pa))/(2*Kk); return Pa/(2*Pa - p) + bv*Pa/(Rv*TT)
Ps = np.logspace(2, 11, 4000)
z400 = Zf(Ps, 1.0e-8, 400); print('b/RT(400)', bv/(Rv*400), 'min Z 400K', z400.min()); assert (z400 > 1).all()
assert sig(bv/(Rv*400), 2) == 1.2e-8
z300 = Zf(Ps, Kv, 300); print('min Z 300K', z300.min(), Ps[z300.argmin()]); assert z300.min() < 1
print('ALL OK')
