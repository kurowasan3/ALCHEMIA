# 検証スクリプト（グループ g11）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-MOL-051/item.json'
import sys, json, re, math
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
eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
assert check_equation('2XeF2 + 2H2O -> 2Xe + 3HF + O2') is not None  # 検査器が誤りを検出できること
def sig(x, n):
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
assert abs(math.sin(math.radians(54.75)) - 0.817) < 5e-4 and abs(math.cos(math.radians(30)) - 0.866) < 5e-4
per = 0.154*0.817; print(per, 17*per); assert sig(per, 3) == 0.126 and sig(17*per, 3) == 2.14
assert abs(math.cos(math.radians((180-109.5)/2)) - math.sin(math.radians(54.75))) < 1e-12
c = 28.4e-3/284/0.100; assert abs(c - 1.00e-3) < 1e-9
N = c*0.250e-3*6.02e23; A = 346/N/1e-14; print(N, A); assert sig(A, 3) == 0.230
t = 15*2*2.50; assert t == 75.0
cosa = 65.0/t; print(cosa, math.degrees(math.acos(cosa))); assert sig(cosa, 3) == 0.867 and abs(math.degrees(math.acos(cosa)) - 30) < 0.5
A2 = 0.200/cosa; print(A2); assert sig(A2, 3) == 0.231
print('OK')
