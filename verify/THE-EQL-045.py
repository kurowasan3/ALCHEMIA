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

it = check_item_equations('THE-EQL-045')
K = 4; Ksp = sp.Rational(18, 10)*sp.Rational(1, 10**10); F = 96500
x = sp.symbols('x')
sols = sp.solve(sp.Rational(2, 100) + x - K*(sp.Rational(1, 10) - x)*(sp.Rational(2, 10) - x), x)
sols = sorted(float(s) for s in sols); print('roots', sols)
xx = [s for s in sols if 0 < s < 0.1]; assert len(xx) == 1; xx = xx[0]
assert sols[1] > 0.1
Q2 = xx*0.1*F; print('Q2', Q2, 'Ag', 0.2-xx); assert r2(Q2, 2) == '2.8e+02' and r2(0.2-xx, 2) == '0.17'
Ag = Ksp/sp.Rational(1, 10); Q = sp.Rational(2, 100)/(sp.Rational(1, 10)*Ag); print('Q', float(Q)); assert r2(Q, 2) == '1.1e+08'
Kp = 1/(K*Ksp); print("K'", float(Kp))
Q4 = 2e-3*F; assert r2(Q4, 2) == '1.9e+02'
Fe3 = K*0.12*(float(Ksp)/0.08); print('Fe3', Fe3); assert r2(Fe3, 2) == '1.1e-09'
# exact solve of the reverse equilibrium to confirm near-completeness
y = sp.symbols('y', positive=True)  # mol/L of Fe3+ reacted
yy = sp.symbols('yy')
roots = sp.Poly((sp.Rational(1, 10) + yy) - Kp*(sp.Rational(2, 100) - yy)*(sp.Rational(1, 10) - yy), yy).nroots(n=40)
s = [r for r in roots if 0 < r < sp.Rational(2, 100)][0]
print('reacted', s, 'left', float(sp.Rational(2, 100) - s))
assert abs(float(sp.Rational(2, 100) - s) - Fe3) / Fe3 < 0.01
print('THE-EQL-045 OK')
