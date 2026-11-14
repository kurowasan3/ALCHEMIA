# 検証スクリプト（東工大 g12）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
ITEM = '/home/claude/chemx/tokodai/out/items/THE-ELC-072/item.json'
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
for key in ('problem','answer','solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>', b)) == 1:
            assert check_equation(b) is None, (b, check_equation(b)); print('式OK', b)
assert it['problem'].count('\\item') == len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])) == 3
F = 96500; Cu = 64; Ag = 108
nCu = R(150,1000)*R(100,1000); nAg = R(200,1000)*R(100,1000)
eCu, eAg = 2*nCu, nAg
print('完全析出の電子', float(eCu), float(eAg)); assert eCu > eAg   # X=CuSO4
# 逆の割り当てなら X<Y で矛盾
assert not (eAg > eCu)
# 実験2
e2 = R(150,100)*965/F; print('e2', float(e2)); assert e2 == R(15,1000)
assert e2 < eCu and e2 < eAg
m2 = e2/2*Cu + e2*Ag; print('質量和', float(m2)); assert abs(float(m2)-2.10) < 1e-9
# (3) 等しい電気量になる体積比と質量比
from sympy import symbols, solve
vx = R(2); vy = R(3)
assert R(300,1000)*vx == R(200,1000)*vy
n = symbols('n', positive=True)
ratio = (n/2*Cu)/(n*Ag); print('Cu:Ag', ratio); assert ratio == R(8,27)
print('ALL OK')
