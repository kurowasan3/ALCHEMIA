# 検証スクリプト INO-AMP-018  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/INO-AMP-018/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(re.findall(r'\((\d+)\)', it['answer'])) == len(nums), '答え中に余分な (n)'
    if n_items: assert nums == [str(i+1) for i in range(n_items)], (n_items, nums)
    else: assert not nums
    print('小問数', n_items, '答え', nums)
check_structure()
check_all_equations()
def bis(f, lo=-16, hi=1):
    """f(h) が h について増加 → f=0 の h（log 空間の二分法）"""
    for _ in range(200):
        m = (lo + hi) / 2
        if f(10**m) > 0: hi = m
        else: lo = m
    return 10**((lo + hi) / 2)
# ---- 個別の検証 ----
need = ['2Al + 3CO2 -> Al2O3 + 3CO', 'CO2 + C -> 2CO', 'C + 2O^2- -> CO2 + 4e-', 'C + O^2- -> CO + 2e-', 'Al^3+ + 3e- -> Al']
eqs = all_equations()
for e in need: assert e in eqs, e
ne = 120.0
co2, co = 30.0*0.800, 30.0*0.200
assert abs(4*co2 + 2*co - 108) < 1e-9
assert ne/4 == 30.0
Al = ne/3 - 2/3*co; print(Al); assert abs(Al-36.0)<1e-9
assert sig(Al*27/1000,3) == 0.972 and sig(Al/(ne/3)*100,3) == 90.0
import sympy as sp
n, y, x = sp.symbols('n y x', positive=True)
ysol = sp.solve(sp.Eq(x, (n/4 - y)/(n/4)), y)[0]
eta = sp.simplify((n/3 - sp.Rational(2,3)*ysol)/(n/3)); print(eta)
assert sp.simplify(eta - (1+x)/2) == 0
# (4)
z = 33.0 - ne/4; yy = 12.0 - 2*z
assert abs(z-3.0)<1e-9 and abs(yy-6.0)<1e-9 and abs(ne/4 - yy - z - 21.0) < 1e-9
Al2 = ne/3 - 2/3*yy; assert abs(Al2-36.0)<1e-9
naive = (1 + 21.0/33.0)/2*100; print(naive); assert sig(naive,3) == 81.8
assert (ne/4 + z)*12 == 396

print('ALL OK')
