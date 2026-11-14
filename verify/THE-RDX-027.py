# 検証スクリプト THE-RDX-027  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-RDX-027/item.json'
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
need = ['2KMnO4 + 3H2SO4 + 5H2O2 -> K2SO4 + 2MnSO4 + 5O2 + 8H2O', '2KMnO4 + 3H2SO4 + 7H2O2 -> K2SO4 + 2MnSO4 + 6O2 + 10H2O', '2H2O2 -> 2H2O + O2', '4KMnO4 + 6H2SO4 + 12H2O2 -> 2K2SO4 + 4MnSO4 + 11O2 + 18H2O']
eqs = all_equations()
for e in need: assert e in eqs, e
import sympy as sp
k = sp.symbols('k')
# 元素ごとの左右差： K, Mn, S, H, O
L = {'K':2, 'Mn':2, 'S':3, 'H':3*2 + (5+2*k)*2, 'O':2*4 + 3*4 + (5+2*k)*2}
R = {'K':2, 'Mn':2, 'S':1+2, 'H':(8+2*k)*2, 'O':4 + 2*4 + (5+k)*2 + (8+2*k)}
for e in L: assert sp.simplify(L[e]-R[e]) == 0, e
# 表の数値
assert (L['H'].subs(k,0), L['O'].subs(k,0), L['H'].subs(k,1), L['O'].subs(k,1)) == (16, 30, 20, 34)
nM = 0.0200*16.0e-3; nO2 = 8.80e-4
kk = 2*nO2/nM - 5; print(kk); assert abs(kk-0.5) < 1e-9
nH = nM*(5+2*kk)/2; c = nH/20.0e-3; print(nH, c); assert sig(c,3) == 4.80e-2
cA = nM*2.5/20.0e-3; assert sig(cA,3) == 4.00e-2
err = (c-cA)/c*100; assert sig(err,3) == 16.7
assert abs(nM*2.5 - 8.00e-4) < 1e-12 and abs(nH - nM*2.5 - 1.60e-4) < 1e-12

print('ALL OK')
