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
from sympy import Rational, nsimplify
it = json.load(open('/home/claude/chem3/out/items/THE-MOL-031/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
assert check_equation('3MnO2 -> Mn3O4 + O2') is None
M, O = 55.0, 16.0
pct = {k: 100*x*M/(x*M+y*O) for k, (x, y) in {'MO': (1,1), 'M3O4': (3,4), 'M2O3': (2,3), 'MO2': (1,2)}.items()}
print(pct)
assert sorted(pct, key=lambda k: -pct[k]) == ['MO', 'M3O4', 'M2O3', 'MO2']
rA = (36.8/O)/(63.2/M); rB = (27.9/O)/(72.1/M)
print('O/M', rA, rB); assert abs(rA-2) < 0.01 and abs(rB-4/3) < 0.01
assert round(pct['MO2'], 1) == 63.2 and round(pct['M3O4'], 1) == 72.1
loss = (3*(M+2*O) - (3*M+4*O)) / (3*(M+2*O)) * 100
print('loss', loss); assert sig(loss, 2) == 12
print('ALL OK')
