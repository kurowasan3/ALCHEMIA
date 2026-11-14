# 検証スクリプト THE-EQL-060  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-EQL-060/item.json'
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
    keys = ['src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig']
    assert all(k in it for k in keys), [k for k in keys if k not in it]
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

from sympy import symbols, solve, Poly, diff, simplify, log
N, F2, P, kA, k1, k2, k3, K1, K2, Fa = symbols('N F2 P kA k1 k2 k3 K1 K2 Fa', positive=True)
Fsol = solve(K1 - P*Fa/(N*F2), Fa)[0]
rates = {'A': kA*N**2*F2, 'B': k1*N*F2, 'C': k2*N*Fsol, 'D': k3*K2*N**2*F2}
def orders(v):
    return tuple(simplify(diff(log(v), x)*x) for x in (N, F2, P))
ok = [m for m, v in rates.items() if orders(v) == (1, 1, 0)]
print({m: orders(v) for m, v in rates.items()}); assert ok == ['B']
assert orders(rates['A']) == orders(rates['D'])
assert orders(rates['C'])[2] == -1
# 2段階の和が全体式になること（原子の数の加算）
from chemlib import parse_formula
from collections import Counter
def side(s):
    c = Counter()
    for t in s.split(' + '):
        t = t.strip(); n = 1
        if t[0].isdigit(): n = int(t[0]); t = t[1:]
        for k, v in parse_formula(t)[0].items(): c[k] += n*v
    return c
def add(*eqs):
    L, Rr = Counter(), Counter()
    for e in eqs:
        a, b = e.split(' -> '); L += side(a); Rr += side(b)
    return L, Rr
for steps in (['NO2 + F2 -> NO2F + F', 'NO2 + F -> NO2F'], ['2NO2 -> N2O4', 'N2O4 + F2 -> 2NO2F']):
    L, Rr = add(*steps)
    # 両辺に共通の中間体を消すと全体式（原子の数）に一致
    assert L - Rr == side('2NO2 + F2') - side('2NO2F') or (L - Rr) == Counter(), (L, Rr)
    tot_L, tot_R = side('2NO2 + F2'), side('2NO2F')
    assert L - tot_L == Rr - tot_R, (L, Rr)
print('2段階の和 = 全体式（中間体は相殺）OK')
print('J133 OK')
