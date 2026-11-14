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

it = check_item_equations('THE-ELC-028')
R = sp.Rational(83, 10)*10**3; T = 300; V = sp.Rational(83, 10)*sp.Rational(1, 1000); F = 96500
# (2)
CdOH2 = sp.Rational(150, 10000) - sp.Rational(200, 10000)/2; Cd = sp.Rational(30, 10000) + sp.Rational(100, 10000)
assert r2(CdOH2, 2) == '0.005' and r2(Cd, 2) == '0.013'
# (3)
q = sp.Rational(4, 1000); O2gen = q/4
n_left = 3000*V/(R*T); print('O2 left', n_left)
r = (O2gen - n_left)/O2gen; print('absorbed', float(r))
assert abs(float(r) - 0.99) < 1e-12
CdOH2_3 = CdOH2 - q/2 + 2*(O2gen - n_left); print('CdOH2 after', float(CdOH2_3))
assert r2(CdOH2_3, 3) == '0.00498'
# (4)
per_e = sp.Rational(1, 2) - 2*r/4
ne = CdOH2_3/per_e; t = ne*F/sp.Rational(2, 10)/3600
print('ne', float(ne), 'hours', float(t)); assert r2(t, 2) == '1.3e+02'
# (5) no reserve: Cd(OH)2 = 0 at overcharge start
absO2 = r*O2gen; e_to_CdOH2 = 2*2*absO2; H2 = (q - e_to_CdOH2)/2
gas = n_left + H2; dp = gas*R*T/V; print('H2', float(H2), 'dp', float(dp))
assert r2(dp, 2) == '9e+03'
# (6) capacities
assert sp.Rational(200, 10000) < 2*Cd
print('THE-ELC-028 OK')
