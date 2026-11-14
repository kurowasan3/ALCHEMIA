# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/THE-ELC-089/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import Rational as R
F=R(965,10)*10**3
n=(R(1014,10000)-R(10,10000))/(2*F); c=n/R(1,10000)
assert sig(float(n),3)==5.20e-7 and sig(float(c),3)==5.20e-3
mg=float(c*R(1,10)*180*1000); assert sig(mg,3)==93.6, mg
asc=float(R(10,10000)/(2*F)/R(1,10000)); assert sig(asc,2)==5.2e-5, asc
cap=float(2*F*R(250,100)*R(1,10**6)); assert abs(cap-0.4825)<1e-9
assert abs((0.4835-0.0010)-cap)<1e-12
print('数値 OK', float(n), float(c), mg, asc, cap)

print('ALL OK')
