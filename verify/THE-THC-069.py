# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/THE-THC-069/item.json'
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
est=6*394+5*286; assert est==3794 and sig(est,3)==3790
qf=est-2820; assert sig(qf,2)==970
dry=R(2820)*1000/162; assert sig(float(dry),3)==17400
water=R(42,10)*80+R(41)*1000/18
net=R(4,10)*dry-R(6,10)*water; assert sig(float(net),2)==5400, float(net)
wl=dry/(dry+water); assert sig(float(wl)*100,2)==87, float(wl)
print('数値OK', est, qf, float(dry), float(water), float(net), float(wl))

print('ALL OK')
