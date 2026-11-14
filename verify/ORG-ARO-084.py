# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/ORG-ARO-084/item.json'
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

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ['src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer', 'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig']:
    assert k in it, k
# ---- 個別の検証 ----

from fractions import Fraction as Fr
H, C, N, O, Na, Cl, Sn = Fr(1), Fr(12), Fr(14), Fr(16), Fr(23), Fr(71, 2), Fr(119)
benz = 6*C + 6*H; hno3 = H + N + 3*O; anil = 6*C + 7*H + N
assert anil == 93
y1 = (benz + hno3 + Fr(3, 2)*Sn + 7*(H + Cl) + (Na + O + H)) / anil
y2 = (benz + hno3 + 3*2*H) / anil
y3 = (benz + 2*Cl + 2*(N + 3*H)) / anil
print(float(y1), float(y2), float(y3))
assert round(float(y1), 1) == 6.6 and round(float(y2), 1) == 1.6 and round(float(y3), 2) == 1.97
assert y2 < y3 < y1
# Cl bookkeeping per 2 aniline: 14 HCl -> 2 in salt + 12 in 3 SnCl4
assert 14 == 2 + 3*4
print('OK')
