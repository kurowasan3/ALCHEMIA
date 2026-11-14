# 検証スクリプト THE-ATM-016  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-ATM-016/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
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
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
import numpy as np
from scipy.integrate import solve_ivp
NA = 6.02e23
NRa = 1.00/226*NA
ratio = 3.82/(1.60e3*365)
NRn = NRa*ratio
V = NRn/NA*22.4e3
assert sig(NRn,2) == 1.7e16 and sig(V,2) == 6.5e-4, (NRn, V)
t99 = 3.82*2/0.301
assert round(t99) == 25
frac = 0.693/(1.60e3*365)*t99
assert sig(frac,2) == 3.0e-5, frac
mRa = 1.00e6/238*1.60e3/4.47e9*226
assert sig(mRa,2) == 0.34, mRa
# 数値積分（日単位）
l1 = np.log(2)/(1.60e3*365); l2 = np.log(2)/3.82
f = lambda t,y: [-l1*y[0], l1*y[0]-l2*y[1]]
s = solve_ivp(f, [0, 60], [NRa, 0], rtol=1e-10, atol=1, dense_output=True)
Nst = s.sol(60)[1]
assert abs(Nst/NRn - 1) < 1e-3
fr = s.sol(25.4)[1]/s.sol(200)[1] if False else s.sol(25.4)[1]/NRn
assert abs(fr-0.99) < 0.001, fr
# (3) 過渡平衡
for TP, TQ, conv in [(10,1,True),(1,10,False)]:
    a, b = np.log(2)/TP, np.log(2)/TQ
    g = lambda t,y: [-a*y[0], a*y[0]-b*y[1]]
    ss = solve_ivp(g, [0, 80], [1, 0], rtol=1e-10, atol=1e-30, dense_output=True)
    r1 = ss.sol(60)[1]/ss.sol(60)[0]; r2 = ss.sol(80)[1]/ss.sol(80)[0]
    if conv:
        assert abs(r1 - 1/9) < 1e-6 and abs(r2-1/9) < 1e-6 and r1 > TQ/TP
    else:
        assert r2 > 100*r1
print('NRn', NRn, 'V mL', V, 't99', t99, 'frac', frac, 'mRa', mRa)

print('ALL OK')
