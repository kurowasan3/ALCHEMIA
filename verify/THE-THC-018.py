# 検証スクリプト THE-THC-018  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-THC-018/item.json'
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
from sympy import Rational as R, Float
check_structure()
eqs = all_equations(); print(eqs)
for e in eqs:
    assert check_equation(e) is None, e
assert '2MgCl -> Mg + MgCl2' in eqs
sub, I1, I2, I3, D, EA, U2 = 148, 738, 1451, 7733, 243, 349, 2526
H2 = sub + I1 + I2 + D - 2*EA - U2; assert H2 == -644
rCl = R(181, 1000)
k = R(U2) * (R(72,1000) + rCl) / (3*2); print('k =', float(k)); assert sig(float(k), 4) == 106.5
U1 = k * 2*1 / (R(100,1000) + rCl); U3 = k * 4*3 / (R(60,1000) + rCl)
print('U1', float(U1), 'U3', float(U3))
assert sig(float(U1), 3) == 758 and sig(float(U3), 3) == 5300
U1r, U3r = 758, 5300
H1 = sub + I1 + R(D, 2) - EA - U1r; H3 = sub + I1 + I2 + I3 + R(3*D, 2) - 3*EA - U3r
print('H1', float(H1), 'H3', float(H3)); assert H1 == R(-199, 2) and sig(float(H3), 3) == 4090
H3raw = sub + I1 + I2 + I3 + R(3*D, 2) - 3*EA - U3; print('H3(丸め前)', float(H3raw))
disp = H2 - 2*H1; print('不均化', float(disp)); assert disp == -445
step1_cost, step1_gain = I2 + R(D,2), EA + (U2 - U1r)
step2_cost, step2_gain = I3 + R(D,2), EA + (U3r - U2)
print(float(step1_cost), step1_gain, float(step2_cost), step2_gain)
assert step1_cost == R(3145,2) and step1_gain == 2117 and step1_cost - step1_gain == H2 - H1 == R(-1089,2)
assert step2_cost == R(15709,2) and step2_gain == 3123 and step2_cost - step2_gain == H3 - H2 == R(9463,2)
print('全検証OK:', ITEM)
