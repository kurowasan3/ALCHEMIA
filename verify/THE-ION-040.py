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
# ---- THE-ION-040 純粋なグリシン水溶液の pH（厳密解と近似式の比較）
from mpmath import mp, mpf, log10, sqrt, findroot
mp.dps = 30
it = load('THE-ION-040'); check_all_equations(it); check_subitems(it)
K1, K2, Kw = mpf('4.0e-3'), mpf('2.5e-10'), mpf('1.0e-14')
def exact(C):
    # 未知数 h, Z(=[G±]) を物質収支と電荷収支から厳密に解く
    def f(lh, lz):
        h, Z = mpf(10)**lh, mpf(10)**lz
        return [h + h*Z/K1 - K2*Z/h - Kw/h, Z*(1 + h/K1 + K2/h) - C]
    sol = findroot(lambda a, b: f(a, b), (mpf(-7), log10(C)))
    return mpf(10)**sol[0]
def formula(C): return sqrt(K1*(K2*C + Kw)/(K1 + C))
res = {}
for C in ['1e-1', '1e-3', '1e-6']:
    C = mpf(C); he, hf = exact(C), formula(C)
    res[str(C)] = -log10(hf)
    print(float(C), 'exact pH', float(-log10(he)), 'formula', float(-log10(hf)))
    assert abs(log10(he) - log10(hf)) < 1e-3
pHs = [round(float(v), 1) for v in res.values()]
assert pHs == [6.0, 6.3, 7.0], pHs
assert abs(-log10(sqrt(K1*K2)) - 6) < 1e-12
# 与えた log 値
for x, v in [(1.03, 0.01), (2.08, 0.32), (4.56, 0.66), (9.62, 0.98)]: assert round(float(log10(x)), 2) == v
# 中間値
assert round(float(K1*(K2*mpf('0.1')+Kw)/(K1+mpf('0.1')))*1e13, 2) == 9.62
assert round(float(K1*(K2*mpf('1e-3')+Kw)/(K1+mpf('1e-3')))*1e13, 2) == 2.08
# (4)
C = mpf('1e-3'); h = formula(C)
Gp, Gm = C*h/K1, K2*C/h
print('G+', float(Gp), 'G-', float(Gm), 'OH-H', float(Kw/h - h))
assert sig(float(Gp), 2) == 1.1e-7 and sig(float(Gm), 2) == 5.5e-7
assert abs((Gp - Gm) - (Kw/h - h)) / abs(Gp - Gm) < 1e-3 and sig(float(Gp - Gm), 3) == -4.34e-7
assert float(Gp/C) < 1e-3 and float(Gm/C) < 1e-3   # [G±]≒C の妥当性
print('THE-ION-040 OK')
