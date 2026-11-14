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
it = json.load(open('/home/claude/chem3/out/items/THE-BND-048/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
data = {'Li': ('bcc', 6.9, 0.351, 4.32e-23), 'Ni': ('fcc', 58.7, 0.352, 4.36e-23), 'Ba': ('bcc', 137.3, 0.502, 1.27e-22), 'Ag': ('fcc', 107.9, 0.409, 6.84e-23)}
r, d = {}, {}
for k, (s, M, a, V) in data.items():
    assert abs((a*1e-7)**3 - V)/V < 0.005, k
    r[k] = (1.73 if s == 'bcc' else 1.41)*a/4
    d[k] = (2 if s == 'bcc' else 4)*M/(6.0e23*V)
print(r); print(d)
assert sorted(r, key=lambda k: -r[k]) == ['Ba', 'Li', 'Ag', 'Ni']
assert sorted(d, key=lambda k: -d[k]) == ['Ag', 'Ni', 'Ba', 'Li']
assert sig(r['Li']/r['Ni'], 2) == 1.2 and sig(d['Li'], 2) == 0.53
print('ALL OK')
