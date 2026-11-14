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
it = json.load(open('/home/claude/chem3/out/items/THE-SOL-049/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
R, T = 8.31e3, 273
n0 = 1.0e5*0.0235/(R*T); n1 = n0*3.0*2.00
print(n0, n1); assert sig(n1, 2) == 6.2e-3
va = n1*R*T/3.0e5; vb = n1*R*T/1.0e5
print(va, vb); assert sig(va, 2) == 0.047 and sig(vb, 2) == 0.14
nN = n0*4.0; nO = 1.0e5*0.0489/(R*T)*1.0
vN = nN*R*T/5.0e5; vO = nO*R*T/5.0e5
print(vN, vO); assert sig(vN, 2) == 0.019 and sig(vO, 2) == 0.0098
vN2 = (n0*8.0)*R*T/1.0e6; print(vN2/vN); assert abs(vN2/vN-1) < 1e-12
print('ALL OK')
