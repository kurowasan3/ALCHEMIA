# 検証スクリプト（グループ g11）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-ELC-092/item.json'
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
# 負極×2＋正極 で全体の反応になることを確認
from collections import Counter
pairs = [('H2 -> 2H+ + 2e-', 'O2 + 4H+ + 4e- -> 2H2O'), ('H2 + 2OH- -> 2H2O + 2e-', 'O2 + 2H2O + 4e- -> 4OH-'),
         ('H2 + O^2- -> H2O + 2e-', 'O2 + 4e- -> 2O^2-'), ('H2 + CO3^2- -> H2O + CO2 + 2e-', 'O2 + 2CO2 + 4e- -> 2CO3^2-')]
def side(s):
    c = Counter()
    for t in s.split(' + '):
        t = t.strip(); m = re.match(r'^(\d*)(.+)$', t); c[m.group(2)] += int(m.group(1) or 1)
    return c
for an, ca in pairs:
    assert check_equation(an) is None and check_equation(ca) is None
    L = Counter(); R = Counter()
    for eq, k in ((an, 2), (ca, 1)):
        l, r = eq.split(' -> ')
        for sp_, n in side(l).items(): L[sp_] += k*n
        for sp_, n in side(r).items(): R[sp_] += k*n
    for sp_ in list(L):
        common = min(L[sp_], R.get(sp_, 0)); L[sp_] -= common; R[sp_] -= common
    L = +L; R = +R; print(dict(L), '->', dict(R))
    assert L == Counter({'H2': 2, 'O2': 1}) and R == Counter({'H2O': 2})
assert check_equation('CO + O^2- -> CO2 + 2e-') is None and check_equation('CO + CO3^2- -> 2CO2 + 2e-') is None
assert check_equation('CO2 + 2OH- -> CO3^2- + H2O') is None
ne = 10.0*3600/9.65e4; nH2 = ne/2; print(ne, nH2, nH2*18.0, nH2*44.0)
assert sig(nH2, 3) == 0.187 and sig(nH2*18.0, 3) == 3.36 and sig(nH2*44.0, 3) == 8.21
print('OK')
