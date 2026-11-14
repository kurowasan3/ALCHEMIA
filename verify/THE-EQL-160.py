# 検証スクリプト（グループ g12）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_all_eqs(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)

def check_parts(it):
    n = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n + 1)], (n, nums)
    print('小問数OK:', n)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

ITEM = '/home/claude/chemk/outh/items/THE-EQL-160/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_all_eqs(it); check_parts(it)
from sympy import symbols, binomial, simplify, Rational as R, nsolve, Symbol
a, c, L0, R0 = symbols('a c L0 R0', positive=True)
Ri = [binomial(4, i) * a**i * R0 for i in range(5)]
Ti = [binomial(4, i) * (c*a)**i * L0 * R0 for i in range(5)]
f = simplify(sum(i*(Ri[i] + Ti[i]) for i in range(5)) / (4 * sum(Ri + Ti)))
fa = (a*(1+a)**3 + L0*c*a*(1+c*a)**3) / ((1+a)**4 + L0*(1+c*a)**4)
xR = simplify(sum(Ri) / sum(Ri + Ti)); xa = (1+a)**4 / ((1+a)**4 + L0*(1+c*a)**4)
assert simplify(f - fa) == 0 and simplify(xR - xa) == 0
print('(1) 式OK')
P = {c: R(1, 100), L0: 1000}
F = lambda av, p=P: float(fa.subs(p).subs(a, av))
f1, f10 = F(1), F(10); print('(2)', f1, f10)
assert (sig(f1, 2), sig(f10, 2)) == (0.017, 0.83)
rel = f10 - f1; rel1 = 10/11 - 1/2; print('放出', rel, rel1, rel/rel1)
assert sig(rel/rel1, 2) == 2.0 and sig(rel, 2) == 0.82
x = Symbol('x', positive=True)
ah = nsolve(((1+x)/(1+x/100))**4 - 1000, x, 3); print('(3) a =', ah); assert sig(float(ah), 2) == 4.9
print('x_R(1), x_R(10) =', float(xa.subs(P).subs(a, 1)), float(xa.subs(P).subs(a, 10)))
assert round(1/float(xa.subs(P).subs(a, 1)) - 1) == 65
f0 = lambda av: float((a*(1+a)**3/((1+a)**4 + L0)).subs(P).subs(a, av))
print('(4)', f0(1), f0(10)); assert (sig(f0(1), 2), sig(f0(10), 2)) == (0.0079, 0.85)
Tsh = lambda av: float((L0*c*a*(1+c*a)**3 / (a*(1+a)**3 + L0*c*a*(1+c*a)**3)).subs(P).subs(a, av))
print('T 構造に結合した割合', Tsh(1), Tsh(10)); assert round(Tsh(1)*100) == 56 and round(Tsh(10)*100) == 1
print('誤差 a=10:', (f0(10) - f10)/f10); assert round((f0(10) - f10)/f10*100) == 2
P2 = {c: R(1, 100), L0: 10000}
ah2 = nsolve(((1+x)/(1+x/100))**4 - 10000, x, 5); print('(5) a =', ah2); assert abs(float(ah2) - 10) < 1e-9
g1, g10 = F(1, P2), F(10, P2); print(g1, g10, g10 - g1)
assert fa.subs(P2).subs(a, 10) == R(1, 2) and sig(g1, 2) == 0.011 and sig(g10 - g1, 2) == 0.49
print('ALL OK')
