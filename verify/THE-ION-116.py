ITEM = '/home/claude/chemk/out/items/THE-ION-116/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
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

f = 5.70
d = {'Cu+':50.0,'CuCl':-119.9,'CuBr':-100.8,'CuI':-69.5,'Cl-':-131.2,'Br-':-104.0,'I-':-51.6,'Ag+':77.1,'CrO4':-727.8,'Ag2CrO4':-641.8}
G = {s: d['Cu+']+d[a]-d[s] for s, a in [('CuCl','Cl-'),('CuBr','Br-'),('CuI','I-')]}
G['Ag2CrO4'] = 2*d['Ag+']+d['CrO4']-d['Ag2CrO4']
lk = {s: -g/f for s, g in G.items()}; print(G, lk)
assert abs(G['CuCl']-38.7) < 1e-9 and [round(lk[s],1) for s in ('CuCl','CuBr','CuI','Ag2CrO4')] == [-6.8,-8.2,-11.9,-12.0]
assert abs(2.303*8.31*298/1000 - 5.70) < 0.01
sCuI = (10**lk['CuI'])**0.5; sAg = (10**lk['Ag2CrO4']/4)**(1/3)
sCuI2 = (10**-11.9)**0.5; sAg2 = (10**-12.0/4)**(1/3)
print(sCuI, sAg, sCuI2, sAg2)
assert sig(sCuI,1) == sig(sCuI2,1) == 1e-6 and sig(sAg,1) == sig(sAg2,1) == 6e-5
s1 = [10**(lk[s]/2) for s in ('CuCl','CuBr','CuI')]; assert s1 == sorted(s1, reverse=True)

print('全検証OK:', ITEM)
