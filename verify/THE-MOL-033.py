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
it = json.load(open('/home/claude/chem3/out/items/THE-MOL-033/item.json', encoding='utf-8'))
eqs = check_all_equations(it); check_parts(it); assert len(eqs) == 5, eqs
H, C, O, Ca = 1.0, 12, 16, 40
f = {'hyd': Ca+2*C+4*O+2*H+O, 'ox': Ca+2*C+4*O, 'carb': Ca+C+3*O, 'oxi': Ca+O}
print(f)
r = [round(100*f[k]/f['hyd'], 1) for k in ('ox', 'carb', 'oxi')]
print(r); assert r == [87.7, 68.5, 38.4]
print('ALL OK')
