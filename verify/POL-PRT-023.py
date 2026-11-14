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
# ---- POL-PRT-023 酸性・塩基性アミノ酸の滴定（厳密な電荷収支で表のデータを再現し，答えを独立に計算）
from mpmath import mp, mpf, log10
mp.dps = 30
it = load('POL-PRT-023'); check_all_equations(it); check_subitems(it)
Kw = mpf('1e-14')
def titr_pH(pks, z0, r, C=mpf('0.1')):
    K = [mpf(10)**(-mpf(p)) for p in pks]
    def g(lh):
        h = mpf(10)**(-lh)
        t = [h**3, h**2*K[0], h*K[0]*K[1], K[0]*K[1]*K[2]]
        ch = sum(ti/sum(t)*(z0 - i) for i, ti in enumerate(t)) * C
        return h + C*mpf(r) + ch - z0*C - Kw/h      # Na+ = rC, Cl- = z0*C
    a, b = mpf('0.5'), mpf('13.5')
    for _ in range(200):
        m = (a + b)/2
        if g(a)*g(m) <= 0: b = m
        else: a = m
    return float(m)
GLU, LYS = (['2.19', '4.25', '9.67'], 1), (['2.18', '8.95', '10.53'], 2)
table = {'X': [3.24, 4.26, 6.96, 9.67], 'Y': [5.58, 8.92, 9.74, 10.56]}
for name, (pks, z0) in [('X', GLU), ('Y', LYS)]:
    sim = [titr_pH(pks, z0, r) for r in [1.0, 1.5, 2.0, 2.5]]
    print(name, [round(x, 3) for x in sim])
    assert all(abs(s - t) <= 0.006 for s, t in zip(sim, table[name])), (sim, table[name])
# 逆に：Y がグルタミン酸型だったら r=1.5 は 4 付近 → 表と矛盾（同定の一意性）
assert abs(titr_pH(GLU[0], 1, 1.5) - 8.92) > 3 and abs(titr_pH(LYS[0], 2, 1.5) - 4.26) > 3
# (1)(2) 表から読む
pK1 = 2.2
pIX = (pK1 + 4.26)/2; pIY = (8.92 + 10.56)/2
print('pI', pIX, pIY); assert round(pIX, 1) == 3.2 and round(pIY, 1) == 9.7
assert abs(table['X'][0] - pIX) < 0.05 and abs(table['Y'][2] - pIY) < 0.05
# (3) pH 4.26 の平均電荷（グルタミン酸，pK = 2.2, 4.26, 9.67）
h = 10**-4.26; K = [10**-2.2, 10**-4.26, 10**-9.67]
t = [h**3, h**2*K[0], h*K[0]*K[1], K[0]*K[1]*K[2]]
avg = sum(ti/sum(t)*(1 - i) for i, ti in enumerate(t))
print('平均電荷', avg); assert round(avg, 2) == -0.49
# pH 6.0 での電荷の符号
def avgq(pH, pks, z0):
    h = 10**-pH; K = [10**-p for p in pks]
    t = [h**3, h**2*K[0], h*K[0]*K[1], K[0]*K[1]*K[2]]
    return sum(ti/sum(t)*(z0 - i) for i, ti in enumerate(t))
assert avgq(6.0, [2.2, 4.26, 9.67], 1) < 0 and avgq(6.0, [2.2, 8.92, 10.56], 2) > 0
# (4)
assert abs((4.26 + 9.67)/2 - 6.965) < 1e-9 and round((2.2 + 8.92)/2, 2) == 5.56
# 近似の確認（(1) 解説の数値）
assert round(10**-4.26 * 1e5, 1) == 5.5 and round(10**(-(14 - 10.56)) * 1e4, 1) == 3.6
for s_ in ['X はグルタミン酸，Y はリシン', '(3) $-0.49$']: in_text(s_, it)
print('POL-PRT-023 OK')
