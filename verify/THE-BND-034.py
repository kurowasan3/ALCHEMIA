# 検証スクリプト THE-BND-034  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-BND-034/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
                eqs.add(b)
    return sorted(eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
from sympy import Rational as R
check_structure()
assert all_equations() == []
K, Ca, Cl = 39, 40, R(71, 2)
comp = {'K': R(900,1000), 'Ca': R(50,1000), 'Cl': 1}
assert comp['K'] + 2*comp['Ca'] == comp['Cl']          # 電気的中性
VN = R(250,1000) * R(1,10**21) * 6 * 10**23; assert VN == 150
# A
A_Cl = 4; A_K = 4*comp['K']; A_Ca = 4*comp['Ca']; A_vac = 4 - A_K - A_Ca
assert (A_K, A_Ca, A_vac) == (R(18,5), R(1,5), R(1,5))
mA = A_K*K + A_Ca*Ca + A_Cl*Cl; print('mA', float(mA))
# B
s = comp['K'] + comp['Ca']; B_K = 4*comp['K']/s; B_Ca = 4*comp['Ca']/s; B_Cl = 4/s
print('B', float(B_K), float(B_Ca), float(B_Cl))
assert [sig(float(x), 3) for x in (B_K, B_Ca, B_Cl)] == [3.79, 0.211, 4.21]
mB = B_K*K + B_Ca*Ca + B_Cl*Cl
# C
C_Cl = 4; C_K = 4*comp['K']; C_Ca = 4*comp['Ca']; C_vac = 4 - C_K
assert C_vac == R(2,5) and (C_K, C_Ca, C_Cl) == (A_K, A_Ca, A_Cl)
mC = C_K*K + C_Ca*Ca + C_Cl*Cl
rho = {k: float(m/VN) for k, m in {'pure': 4*(K+Cl), 'A': mA, 'B': mB, 'C': mC}.items()}
print(rho)
assert {k: sig(v, 3) for k, v in rho.items()} == {'pure': 1.99, 'A': 1.94, 'B': 2.04, 'C': 1.94}
# 実際の a^3 でも結論が同じか
a3 = 0.630**3 * 1e-21 * 6.0e23
print({k: round(float(m)/a3, 3) for k, m in {'pure': 4*(K+Cl), 'A': mA, 'B': mB}.items()})
meas = 1.94
assert abs(rho['A'] - meas) < 0.01 and abs(rho['B'] - meas) > 0.05
hole = 1.732*0.630/4 - 0.181; print('hole', hole); assert sig(hole, 2) == 0.092 and hole < 0.100
print('全検証OK:', ITEM)
