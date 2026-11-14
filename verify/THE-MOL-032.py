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
from sympy import symbols, solve, Eq
it = json.load(open('/home/claude/chem3/out/items/THE-MOL-032/item.json', encoding='utf-8'))
eqs = check_all_equations(it); check_parts(it); assert len(eqs) == 2
x, y, a, b = symbols('x y a b')
print(solve([Eq(20*x, 2*10), Eq(20*y, 2*30)], [x, y]))
print(solve([Eq(20*x, 2*20), Eq(20*y, 2*40)], [x, y]))
s = solve([Eq(a+b, 30), Eq(2*a+3*b, 80)], [a, b]); print(s)
assert s[a] == 10 and s[b] == 20
assert s[a]/2 + s[b] == 25 and 3*s[a]/2 + 2*s[b] == 55
cands = [(i, 4-i) for i in range(1, 4)]; print('x+y=4 candidates', cands); assert len(cands) > 1
print('ALL OK')
