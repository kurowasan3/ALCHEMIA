# 検証スクリプト（グループ g08）  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-ION-086/item.json'
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if should_check(b):
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n_eq += 1
print('反応式のつり合い OK：', n_eq, '式')
# 小問数と答えの数の一致
m_en = re.search(r'\\begin\{enumerate\}(.*)\\end\{enumerate\}', it['problem'], re.S)
nq = m_en.group(1).count('\\item') if m_en else 0
assert '\\begin{itemize}' not in it['problem']
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問数', nq, '答えの数', na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


from sympy import symbols, solve, Eq, nsimplify
import math
Ka = 5.6e-10; H = 0.60; R = 8.31; T = 300
pKa = -math.log10(Ka); assert abs(pKa-9.25) < 0.003
al = lambda pH: Ka/(Ka+10**(-pH))
a = [al(pH) for pH in (7.25, 9.25, 11.25)]
# 問題文の近似（pKa=9.25 ちょうど）
a2 = [1/101, 0.5, 100/101]
print(a, a2); assert [sig(x,2) for x in a2] == [0.0099, 0.5, 0.99]
p, c0 = symbols('p c0', positive=True)
def f_of(alpha, V):
    sol = solve(Eq(H*p/alpha + p*V/(R*T), c0), p)[0]
    return (sol*V/(R*T)/c0).simplify()
f1 = float(f_of(100/101, 1000)); f2 = float(f_of(1/101, 1000))
print(f1, f2)
assert sig(f1,2) == 0.40 and sig(f2,2) == 0.0066
n = 1
while (1-f1)**n > 0.10: n += 1
print('回数', n, (1-f1)**4, (1-f1)**5); assert n == 5
V1 = 9*(H/(100/101))*R*T; print('一回で必要', V1); assert sig(V1,2) == 14000
assert V1 > 5*1000

print('ALL OK')
