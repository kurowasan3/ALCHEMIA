# 検証スクリプト THE-EQL-087（グループ g09）  python3 -I この.py
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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    print('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem6/outh/items/THE-EQL-087/item.json')
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
# 種: E*（標識エステル）, E0（非標識エステル）, T*（18OH と 16OH を 1 個ずつ）, T0（標識なし）
def run(p, k=1.0, kT=20000.0, tmax=3.0):
    kr, kf = p*kT, (1-p)*kT
    def f(t, z):
        Es, E0, Ts, T0 = z
        return [-k*Es + kr*Ts/2,
                -k*E0 + kr*Ts/2 + kr*T0,
                k*Es - (kr+kf)*Ts,
                k*E0 - (kr+kf)*T0]
    ts = np.linspace(0.3, tmax, 40)
    sol = solve_ivp(f, (0, tmax), [1, 0, 0, 0], t_eval=ts, method='LSODA', rtol=1e-10, atol=1e-13)
    Es, E0 = sol.y[0], sol.y[1]
    y = Es + E0; x = Es / y
    return np.log(x) / np.log(y)
for p in (0.2, 0.392, 0.7):
    r = run(p)
    print(p, r.min(), r.max(), 'theory', p/(2*(1-p)))
    assert abs(r.mean() - p/(2*(1-p))) < 3e-3 and r.max()-r.min() < 3e-3
P = sp.symbols('p', positive=True)
ratio = sp.log(sp.Rational(8,10))/sp.log(sp.Rational(1,2))
print('比', float(ratio), ' 与えた値で', -0.223/-0.693)
psol = sp.solve(sp.Eq(P/(2*(1-P)), ratio), P)[0]
print('p =', float(psol)); assert round(float(psol), 2) == 0.39
x25 = float(sp.Rational(1,4)**ratio); print('x(y=0.25)=', x25); assert round(x25, 2) == 0.64
print('e^-0.446 =', np.exp(-0.446), ' 0.3218*ln0.25=', 0.3218*np.log(0.25))
print('OK')
