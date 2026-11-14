# 検証スクリプト（グループ g08）  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/THE-GAS-048/item.json'
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


# ファンデルワールス式で低温・中圧の CO2 の圧縮因子を確認
from sympy import symbols, nsolve
a, b, R = 0.364, 4.27e-5, 8.314   # Pa m6/mol2, m3/mol
P, T = 1.0e6, 300.0
x = symbols('x')
v = nsolve((P - (R*T/(x*1e-3-b) - a/(x*1e-3)**2))/P, x, R*T/P*1e3)*1e-3
Z = P*v/(R*T); print('Z =', Z); assert Z < 1
ans = it['answer']
assert '(a)，(b)，(d)，(e)' in ans

print('ALL OK')
