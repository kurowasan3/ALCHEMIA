ITEM = '/home/claude/chem2/out/items/THE-ION-046/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
assert check_equation('Ag2CrO4(s) <=> 2Ag+ + CrO4^2-') is None
from sympy import Rational as R, symbols, nsolve, solve, real_roots
s = symbols('s', positive=True)
K = R(4, 10**12)
s1 = [r for r in solve(4*s**3 - K, s)][0]; print('(1)', float(s1)); assert sig(float(s1), 2) == 1.0e-4
s2 = float(nsolve((2*s)**2 * (R(1,100) + s) - K, s, 1e-5)); print('(2)', s2); assert sig(s2, 2) == 1.0e-5
s3 = float(nsolve((R(1,100) + 2*s)**2 * s - K, s, 4e-8)); print('(3)', s3); assert sig(s3, 2) == 4.0e-8
assert s2 / 0.010 < 0.0011 and 2*s3 / 0.010 < 1e-5

print('全検証OK:', ITEM)
