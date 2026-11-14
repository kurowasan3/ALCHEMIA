ITEM = '/home/claude/chem2/out/items/ORG-STR-016/item.json'
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
eqs = check_all_equations(); assert 'CuO + CO -> Cu + CO2' in eqs
from sympy import Rational as R
def analyze(m, a, b):
    C = b * R(12, 44); H = a * R(2, 18); O = m - C - H
    return C, H, O
C, H, O = analyze(R(22), R(18), R(44))
assert (C, H, O) == (12, 2, 8)
mol = [C/12, H/1, O/16]; r = [v / min(mol) for v in mol]; assert r == [2, 4, 1]
assert 2*(2*12 + 4*1 + 16) == 88
C2, H2, O2 = analyze(R(22), R(18), R(33))
assert C2/22 < C/22 and H2 == H and O2/22 > O/22
water = R('23.8') - 12 * R(44, 24); print('(3) water', water); assert water == R('1.8')
assert R(18) + water == R('19.8')
C3, H3, O3 = analyze(R('23.8'), R('19.8'), R(44)); print(C3, H3, O3)
assert H3 > H and O3 > O and C3/R('23.8') < C/22
C4, H4, O4 = analyze(R(22), R(19), R(45)); assert C4 > C and H4 > H

print('全検証OK:', ITEM)
