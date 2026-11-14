# 検証スクリプト THE-EQL-074  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/THE-EQL-074/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, formula_mass
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
from sympy import symbols, solve, Rational as R
x = symbols('x', positive=True)
def pco(S, K):
    return [s for s in solve(x**2 - K*(S - x), x) if s > 0][0]
n1 = 1 + 1/R(60, 100); print('(1)', float(n1)); assert sig(float(n1), 2) == 2.7
c2 = pco(R(40, 100), R(90, 100)); print('(2)', c2, R(40,100)-c2)
assert c2 == R(3, 10) and (R(40,100)-c2)/c2 < R(60,100)
c3 = pco(1, R(90, 100)); print('(3)', c3, 1-c3)
assert c3 == R(3, 5) and (1-c3)/c3 > R(60,100) and sig(float((1-c3)/c3), 2) == 0.67
c4 = pco(R(40, 100), R(8, 100)); r4 = (R(40,100)-c4)/c4; print('(4)', float(c4), float(r4))
assert sig(float(r4), 2) == 1.8 and r4 > R(70, 100)
# 実際の温度依存式との整合（概数の根拠）
import math
for T in (973.15, 873.15):
    K1 = 10**(688/T - 0.90); Kp2 = 10**(-8916/T + 9.11) * 1.013e5
    print(T, round(K1, 3), '%.3g' % Kp2)
assert 'FeO(s) + CO(g) <=> Fe(s) + CO2(g)' in all_equations() and '2CO -> C + CO2' in all_equations()
assert len(all_equations()) == 3

print('ALL OK')
