# 検証スクリプト THE-SOL-043  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-SOL-043/item.json'
import sys, json, re
from math import log10, floor
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
def check_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '->', 'OK' if r is None else r)
        assert r is None, (e, r)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def close(a, b, rel=1e-9):
    assert abs(a - b) <= rel * max(abs(a), abs(b)), (a, b)
check_structure()
from sympy import symbols, solve, Eq, nsolve
print('反応式', check_equations())
R, T = 8.3e3, 300; RT = R*T
# 作問に使った真の値で測定値を再現
k_true, n_true = 3.0e-7, 0.100
P = lambda V: n_true/(V/RT + k_true)
print('再現 P(2.00)', P(2.0), 'P(0.500)', P(0.5)); assert sig(P(2.0), 3) == 9.06e4 and sig(P(0.5), 3) == 2.00e5
# (1) 測定値から逆算
P1, V1, P2, V2 = 9.06e4, 2.00, 2.00e5, 0.500
print('PV', P1*V1, P2*V2, 'ボイルなら', P1*V1/V2); assert P2*V2 < P1*V1
k, n = symbols('k n')
s = solve([Eq(n, P1*(V1/RT + k)), Eq(n, P2*(V2/RT + k))], [k, n])
kk, nn = float(s[k]), float(s[n]); print('(1) k', kk, 'n', nn); assert sig(kk, 2) == 3.0e-7 and sig(nn, 2) == 0.10
# (2) 2 桁の値を使用
k2, n2 = 3.0e-7, 0.10
plim = n2/k2; print('(2) 極限', plim); assert sig(plim, 2) == 3.3e5
Vg = (n2/3.0e5 - k2)*RT; print('(2) V', Vg, '気体', 3.0e5*Vg/RT, '溶解', k2*3.0e5); assert sig(Vg, 2) == 0.083
# (3)
nN, nC = symbols('nN nC')
eqs = [Eq(nN*RT/V + nC/(V/RT + k2), Pm) for V, Pm in ((2.00, 1.076e5), (0.500, 3.488e5))]
s3 = solve(eqs, [nN, nC]); a, b = float(s3[nN]), float(s3[nC]); print('(3) N2', a, 'CO2', b)
assert sig(a, 2) == 0.050 and sig(b, 2) == 0.050
# 測定値の再現（真値 0.0500, 0.0500）
for V in (2.0, 0.5):
    print('  再現', V, 0.05*RT/V + 0.05/(V/RT + k2))
pN, pC = a*RT/0.5, b/(0.5/RT + k2); frac = pC/(pN + pC); print('(3) 気相 CO2 モル分率', frac); assert sig(frac, 2) == 0.29
print('すべて一致')
