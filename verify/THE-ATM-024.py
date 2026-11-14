# 検証スクリプト THE-ATM-024  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-ATM-024/item.json'
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
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', b): continue  # 一般式（validate と同じ除外）
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式', n, '本 OK')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure(); check_all_equations()
import sympy as sp
p, q = sp.symbols('p q')
sols = sp.solve([p*q-sp.Rational(30,100), (1-p)*(1-q)-sp.Rational(15,100)], [p,q])
print(sols); assert set(sols) == {(sp.Rational(2,5), sp.Rational(3,4)), (sp.Rational(3,4), sp.Rational(2,5))}
P, Q = sp.Rational(2,5), sp.Rational(3,4)
assert P*(1-Q)+(1-P)*Q == sp.Rational(55,100)
nS = 0.0500*10.0e-3; nX = sp.symbols('nX')
NX = float(sp.solve(sp.Rational(1,2)*nX/(nX+nS) - P, nX)[0]); assert sig(NX,3) == 2.00e-3
pp = 0.5*NX/(NX+nS); dist = [pp*0.75, pp*0.25+(1-pp)*0.75, (1-pp)*0.25]
assert all(abs(a-b) < 1e-9 for a,b in zip(dist, [0.30,0.55,0.15]))
nat = [0.5*0.75, 0.5*0.25+0.5*0.75, 0.5*0.25]; assert nat == [0.375, 0.5, 0.125]
print('ALL OK')
