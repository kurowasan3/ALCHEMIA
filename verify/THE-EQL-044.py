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

it = check_item_equations('THE-EQL-044')
y, n = sp.symbols('y n', positive=True)
C = n; O = y*n + 2*(1-y)*n; e = 2*O; Al = e/3
ratio = sp.simplify(C/Al); print('C/Al =', ratio)
assert sp.simplify(ratio - 3/(2*(2-y))) == 0
f = lambda yy: sp.Integer(3)/(2*(2-sp.sympify(yy)))
assert f(0) == sp.Rational(3, 4) and f(1) == sp.Rational(3, 2)
def solve_y(Kp, P):
    return [s for s in sp.solve(y**2*P/(1-y) - Kp, y) if 0 < s < 1][0]
y2 = solve_y(5*10**4, 10**5); print('y2', y2); assert y2 == sp.Rational(1, 2)
C2 = f(y2)*54000/27*12/1000; print('C kg', C2); assert r2(C2, 2) == '24'
y3 = solve_y(5*10**4, 2*10**5); print('y3', float(y3), 'C/Al', float(f(y3)))
assert r2(y3, 2) == '0.39' and r2(f(y3), 2) == '0.93' and f(y3) < f(y2)
rat = sp.Rational(370, 12)/sp.Rational(1000, 27)
y4 = sp.solve(f(y) - rat, y)[0]; Kp4 = y4**2*10**5/(1-y4)
print('y4', float(y4), 'Kp4', float(Kp4)); assert r2(Kp4, 2) == '4.9e+03'
print('THE-EQL-044 OK')
