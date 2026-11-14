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
it = json.load(open('/home/claude/chem3/out/items/THE-GAS-036/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
q, n, M = 500, 0.250, 128
m = n*M
cs = q*4.0/(m*60); cl = q*4.0/(m*40); H = q*10/n/1000
print(cs, cl, H)
assert sig(cs, 2) == 1.0 and sig(cl, 2) == 1.6 and sig(H, 2) == 20
# 0.500 mol
t_start = (2*m*cs*60)/q; t_melt = 2*n*H*1000/q
print(t_start, t_start + t_melt); assert abs(t_start-8.0) < 1e-9 and abs(t_start+t_melt-28) < 1e-9
import xml.etree.ElementTree as ET
ET.parse('/home/claude/chem3/out/items/THE-GAS-036/problem.svg')
print('ALL OK')
