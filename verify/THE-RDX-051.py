# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/THE-RDX-051/item.json'
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

h=6.63e-34; c=3.00e8; NA=6.02e23; E=4.4e-19
lam=h*c/E; assert sig(lam*1e9,2)==450
e500=h*c/500e-9; assert e500<E and sig(e500,2)==4.0e-19
e450=h*c/450e-9; assert e450>=E and sig(e450,3)==4.42e-19
n=0.50*3600/e450/NA; assert sig(n,2)==6.8e-3
frac=2*1.0e-4/n*100; assert sig(frac,2)==3.0
o2=2*1.0e-4/4; assert abs(o2-5.0e-5)<1e-12
print('数値OK', lam, e500, e450, n, frac)

print('ALL OK')
