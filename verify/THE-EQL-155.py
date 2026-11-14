# 検証スクリプト THE-EQL-155  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-EQL-155/item.json'
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|\$', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式の検査', n, '本 OK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

from sympy import symbols, Rational as Q, solve, Eq
free = [30, 38, 12, 20]; obs = [40, 40, 6, 14]
phi, a, b, c, dd = symbols('phi a b c dd')
cav = [a, b, c, dd]
eqs = [Eq(obs[i], (1 - phi)*free[i] + phi*cav[i]) for i in range(4)] + [Eq(a + b + c + dd, 100)]
sol = solve(eqs + [Eq(c, 0)], [phi, a, b, c, dd], dict=True)[0]
print(sol); assert sol[phi] == Q(1, 2) and [sol[a], sol[b], sol[dd]] == [50, 42, 8]
assert Q(1, 2)*(1 + Q(30, 100)) == Q(65, 100)
ex = Q(30, 100)*40; ee_c = ex / (sol[phi]*sol[a]); assert ee_c == Q(12, 25)
phi2 = 1 - Q(3, 12); assert phi2 == Q(3, 4)
pred = [(1 - phi2)*free[i] + phi2*[50, 42, 0, 8][i] for i in range(4)]
print(pred); assert pred == [45, 41, 3, 11]
ee2 = phi2*50*ee_c / pred[0]; print(ee2); assert ee2 == Q(2, 5)
# (6) R=0 なしでは不定
gen = solve(eqs, [phi, a, b, c, dd], dict=True)
print('without R=0', gen)
assert len(gen) == 1 and any(len(v.free_symbols) > 0 for v in gen[0].values())

print('ALL OK')
