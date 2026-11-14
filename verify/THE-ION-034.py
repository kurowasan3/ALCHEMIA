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

it = check_item_equations('THE-ION-034')
Ksp = sp.Rational(1, 10**33); K = 30; h = sp.symbols('h', positive=True)
S = Ksp/h**3 + K*h
hmin = sp.solve(sp.diff(S, h), h)[0]; print('hmin', hmin, float(hmin))
pHmin = 14 + sp.log(hmin, 10); assert round(float(pHmin), 1) == 5.5
Smin = S.subs(h, hmin); print('Smin', float(Smin)); assert r2(Smin, 2) == '1.3e-07'
assert sp.simplify((Ksp/hmin**3)/(K*hmin) - sp.Rational(1, 3)) == 0
# (3) exact boundaries
c = sp.Rational(1, 100)
import mpmath as mp
g = lambda pOH: float(Ksp)/10**(-3*pOH) + K*10**(-pOH) - 0.01
roots = [10**(-mp.findroot(g, (8.0, 12.0), solver='bisect')), 10**(-mp.findroot(g, (2.0, 8.5), solver='bisect'))]
roots = sorted(float(r) for r in roots); print('roots', roots)
pHs = [14 + float(sp.log(r, 10)) for r in roots]; print('pH bounds', pHs)
assert [round(p, 1) for p in pHs] == [3.7, 10.5]
# (4)
Mg = sp.Rational(18, 10)*sp.Rational(1, 10**11)
hMg = sp.sqrt(Mg/sp.Rational(1, 10**5)); pHMg = 14 + sp.log(hMg, 10); print('pH Mg', float(pHMg))
assert round(float(pHMg), 1) == 11.1 and float(pHMg) > pHs[1]
print('S at pH Mg', float(S.subs(h, hMg)))
print('THE-ION-034 OK')
