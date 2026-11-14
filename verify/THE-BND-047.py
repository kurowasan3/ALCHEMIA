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

def check_all_equations(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts(it):
    np_ = it['problem'].count('\\item')
    na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), (np_, na)
    print('小問数OK:', np_)
import math
it = json.load(open('/home/claude/chem3/out/items/THE-BND-047/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
a, rCl = 0.411, 0.181
for s3 in (1.73, math.sqrt(3)):
    r = s3*a/2 - rCl; print('rCs', r); assert round(r, 3) == 0.175 or round(r+1e-9, 3) == 0.175
assert 2*rCl < a
d = (132.9+35.5)/(6.0e23*0.0694e-21); print('d', d); assert sig(d, 2) == 4.0
assert abs(0.411**3 - 0.0694) < 1e-4
print('ALL OK')
