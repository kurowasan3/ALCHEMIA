# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/ORG-ARO-063/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
print('小問の数', np_)

from math import log10, floor, sqrt
def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import symbols, diff, solve, sqrt as ssqrt, simplify, Rational
h, Ka, Kx = symbols('h K_a K_x', positive=True)
g = h*Ka/((h+Kx)*(h+Ka))
crit = solve(diff(g, h), h); print('極値', crit)
assert len(crit) == 1 and simplify(crit[0] - ssqrt(Ka*Kx)) == 0
gmax = simplify(g.subs(h, ssqrt(Ka*Kx))); print('最大値', gmax)
# 速さの等しい2点：h1*h2 = Ka*Kx（一意）
h1, h2 = symbols('h1 h2', positive=True)
eq = simplify((g.subs(h, h1) - g.subs(h, h2)) * (h1+Kx)*(h1+Ka)*(h2+Kx)*(h2+Ka) / (Ka))
print('等速条件', eq.factor())
KA = 1.0e-10
KX = 10**-8.0 * 10**-11.0 / KA; print('Kx', KX); assert abs(KX/1e-9 - 1) < 1e-9
gn = lambda pH, K: (10**-pH)*K/((10**-pH + KX)*(10**-pH + K))
assert abs(gn(8.0, KA) - gn(11.0, KA)) < 1e-15
pHopt = -log10(sqrt(KA*KX)); print('pHopt', pHopt); assert round(pHopt, 1) == 9.5
gm = gn(pHopt, KA); print('最大割合', gm); assert sig(gm*100, 2) == 5.8
assert abs(gm - 1/(1+3.16)**2) < 2e-4
r13 = gn(13.0, KA)/gm; print('pH13/最大', r13); assert sig(r13, 2) == 1.7e-3
# 数値で全域を走査して最大を確認
best = max((gn(p/100, KA), p/100) for p in range(400, 1400)); print('走査最大', best); assert abs(best[1] - 9.5) < 0.011
KB = 1.0e-5
pHa = -log10(sqrt(KB*KX)); ga = gn(pHa, KB); print('アミン', pHa, ga)
assert round(pHa, 1) == 7.0 and sig(ga*100, 2) == 98
best = max((gn(p/100, KB), p/100) for p in range(400, 1400)); assert abs(best[1] - 7.0) < 0.011
print('各成分の割合（フェノール，最適）', 1/(1+sqrt(KX/KA)))
assert round(1/(1+sqrt(KX/KA)), 2) == 0.24
print('ORG-ARO-063 OK')
