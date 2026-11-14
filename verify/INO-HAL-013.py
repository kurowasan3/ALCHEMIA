import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp

def ce_blocks(s):
    out = []; i = 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k = j + 4; d = 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k

def check_item_equations(fid):
    it = json.load(open(f'/home/claude/chem/out2/items/{fid}/item.json', encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for e in ce_blocks(it[key]):
            arrows = len(re.findall(r'->|<=>|<-', e))
            if arrows != 1 or '\\boxed' in e: continue
            lhs, rhs = re.split(r'<=>|->', e)
            if not lhs.strip() or not rhs.strip(): continue
            r = check_equation(e)
            assert r is None, (key, e, r)
            n += 1
    print(f'{fid}: {n} equations balanced')
    return it

def r2(x, n):  # 有効数字 n 桁の文字列
    return f'{float(x):.{n}g}'

it = check_item_equations('INO-HAL-013')
kH = 6.0e-7; K1 = 4.5e-4; Ka = 3.0e-8; Kw = 1e-14
def solve_pure(p):
    c = kH*p
    H, Cl, HClO, ClO, OH = sp.symbols('H Cl HClO ClO OH', positive=True)
    eqs = [H*Cl*HClO - K1*c, H*ClO - Ka*HClO, H*OH - Kw, H - Cl - ClO - OH, Cl - HClO - ClO]
    x0 = (K1*c)**(1/3)
    s = sp.nsolve(eqs, [H, Cl, HClO, ClO, OH], [x0, x0, x0, x0*1e-6, Kw/x0])
    return c, [float(v) for v in s]
c1, s1 = solve_pure(1.0e5); print(c1, s1)
assert r2(c1, 2) == '0.06' and r2(s1[2], 2) == '0.03'
pH1 = -float(sp.log(s1[0], 10)); print('pH', pH1); assert round(pH1, 1) == 1.5
tot1 = c1 + s1[2] + s1[3]; assert r2(tot1, 2) == '0.09'
c2, s2 = solve_pure(1.25e4); tot2 = c2 + s2[2] + s2[3]; print('ratio', tot2/tot1); assert r2(tot2/tot1, 2) == '0.25'
print('fractions', s1[2]/tot1, s2[2]/tot2)
def bleach(pH):
    H = 10**(-pH); Cl = 0.10
    a = H*Cl/K1; b = Ka/H
    HClO = 1e-2/(1 + a + b); return a*HClO, a*HClO/kH
c7, p7 = bleach(7.0); c2_, p2 = bleach(2.0); print(c7, p7, c2_, p2, p2/p7)
assert r2(p7, 2) == '0.28' and r2(c2_, 2) == '0.0069' and r2(p2, 2) == '1.1e+04'
naive = 1e-2*1e-2*0.10/K1; print('naive', naive); assert naive > 1e-2
print('INO-HAL-013 OK')
