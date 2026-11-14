import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def load(item_id):
    return json.load(open(f'/home/claude/chem/out2/items/{item_id}/item.json', encoding='utf-8'))

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
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return len(eqs)

def check_subitems(it):
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    if n_items:
        assert [int(x) for x in nums] == list(range(1, n_items + 1)), (n_items, nums)
    print('小問数OK:', n_items)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def in_text(val_str, it, key='answer'):
    assert val_str in it[key], (val_str, key)
# ---- THE-ION-041 グリシン緩衝液（厳密な電荷収支による数値解と近似の比較）
from mpmath import mp, mpf, log10, sqrt
mp.dps = 30
it = load('THE-ION-041'); n = check_all_equations(it); check_subitems(it)
K1, K2, Kw = mpf('4.5e-3'), mpf('2.5e-10'), mpf('1e-14')
assert round(float(-log10(K1)), 2) == 2.35 and round(float(-log10(K2)), 2) == 9.60
def exact(C, Cl, Na):
    # G+, G±, G-, H+, OH- すべてを含む厳密解（二分法）
    def g(lh):
        h = mpf(10)**(-lh); s = h*h + h*K1 + K1*K2
        return h + Na + C*h*h/s - Cl - C*K1*K2/s - Kw/h
    a, b = mpf(0), mpf(14)
    for _ in range(200):
        m = (a + b)/2
        if g(a)*g(m) <= 0: b = m
        else: a = m
    return m
C = mpf('0.100')
# (1)
naive = -log10(K1) + log10(mpf('0.080')/mpf('0.020'))
ex1 = exact(C, mpf('0.020'), 0)
print('(1)', float(naive), float(ex1))
assert round(2.35 + 2*0.30, 2) == 2.95 and round(float(ex1), 2) == 2.98
# 簡略式（G-, OH- 無視）の二次方程式の根 1.05e-3
h1 = (-mpf('0.0845') + sqrt(mpf('0.0845')**2 + 4*mpf('9.0e-5')))/2
assert round(float(h1)*1e3, 2) == 1.05
# (2)
c_need = C/2 + K1; V = c_need*100/mpf('0.20')
print('(2) V', float(V), 'pH(V=25)', float(exact(C, mpf('0.050'), 0)))
assert abs(float(V) - 27.25) < 1e-9 and sig(float(V), 2) == 27
assert abs(float(exact(C, c_need, 0)) - float(-log10(K1))) < 1e-4
assert round(float(exact(C, mpf('0.050'), 0)), 2) == 2.41
# (3)
d1 = exact(mpf('1e-3'), mpf('5e-4'), 0); d2 = exact(mpf('1e-3'), 0, mpf('5e-4')); u2 = exact(C, 0, mpf('0.050'))
print('(3)', float(d1), float(d2), 'undiluted NaOH side', float(u2))
assert round(float(d1), 2) == 3.38 and round(float(d2), 2) == 9.54 and round(float(u2), 2) == 9.60
assert abs(float(d1) - 2.41) > abs(float(d2) - float(u2))
assert round(float(d1) - 2.41, 2) == 0.97
h3 = mpf(10)**(-d1); print('G+', float(mpf('5e-4') - h3), 'G±', float(mpf('5e-4') + h3))
assert sig(float(mpf('5e-4') - h3), 2) == 8.5e-5 and sig(float(mpf('5e-4') + h3), 2) == 9.2e-4
# 与えた数値表
for x, v in [(2, 0.30), (1.05, 0.021), (1.15, 0.061), (3.86, 0.587), (4.15, 0.618)]:
    assert abs(float(log10(x)) - v) < 0.0015, (x, float(log10(x)))
for x, v in [(75, 8.66), (38.7, 6.22), (34, 5.83)]: assert round(float(sqrt(x)), 2) == v
# (4)
c0 = 0.86*4.5e-3/0.01; print('(4) approx', c0)
# 厳密：pH - pK1 = 0.01 となる c0（h を厳密に解く）
def err(c0):
    c0 = mpf(c0); # G+ 仕込み c0, G± 仕込み c0 → C = 2c0, Cl = c0
    return exact(2*c0, c0, 0) + log10(K1)
lo, hi = mpf('0.05'), mpf('5')
for _ in range(100):
    mid = (lo + hi)/2
    if err(mid) > mpf('0.01'): lo = mid
    else: hi = mid
print('(4) exact c0 for error 0.01:', float(mid))
assert sig(c0, 1) == 0.4 and sig(float(mid), 1) == 0.4 and round(float(mid), 2) == 0.38
for s_ in ['(1) (i) 2.95　(ii) 2.98', '$V=27\\,\\mathrm{mL}$', '(3) 3.38', '(4) $0.4\\,\\mathrm{mol/L}$']: in_text(s_, it)
print('THE-ION-041 OK')
