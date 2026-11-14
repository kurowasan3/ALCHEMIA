# 検証スクリプト（グループ g08）  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/THE-ION-081/item.json'
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


from sympy import symbols, nsolve, log, Rational
import math
Kb, Kw, C = 1.8e-5, 1.0e-14, 0.050
Ka = Kw/Kb
def ph_exact(nNH3, nHCl=0.050, extra_acid=0.0):
    # 電荷: [H+] + [NH4+] = [Cl-] + [OH-]
    Cl = nHCl + extra_acid; N = nNH3
    f = lambda h: h + N*h/(h+Ka) - Cl - Kw/h
    lo, hi = 1e-14, 1.0
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    return -math.log10(mid)
p = [ph_exact(0.020), ph_exact(0.050), ph_exact(0.100)]
print(p)
assert [round(x,1) for x in p] == [1.5, 5.3, 9.3]
# 近似
assert round(2-math.log10(3),1) == 1.5
assert round(0.5*(14-4.74+1.30),2) == 5.28
# 緩衝能: 0.001 mol の強酸を加えたときの pH 変化
dp = [abs(ph_exact(n, extra_acid=0.001)-ph_exact(n)) for n in (0.020,0.050,0.100)]
print('pH変化', dp); assert dp[2] < 0.03 and dp[1] > 2.0
print('(2)+HCl後 pH', ph_exact(0.050, extra_acid=0.001))

print('ALL OK')
