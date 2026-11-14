# 検証スクリプト INO-ALK-019  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/INO-ALK-019/item.json'
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
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
from sympy import Rational as R, Min, symbols, solve, Piecewise
check_structure(); check_all_equations()
assert len(all_equations()) == 4, all_equations()  # 異なる式 4 種（Ca(HCO3)2 の分解は重複）
nCa = R(200,10000)*R(1,2); assert nCa == R(1,100)
def m(n):  # 沈殿の質量 g
    if n <= nCa: return 100*n
    return 100*max(R(0), 2*nCa - n)
assert m(nCa) == 1 and m(2*nCa) == 0
assert 2*nCa*R(224,10) == R(448,1000)
sols = [n for n in [R(k,100000) for k in range(0, 2501)] if m(n) == R(7,10)]
print('0.700 g の解', [float(s) for s in sols]); assert sols == [R(7,1000), R(13,1000)]
# (4) 0.0250 mol 吸収→Ca はすべて Ca(HCO3)2 → 加熱で CaCO3 0.0100 mol
assert R(25,1000) >= 2*nCa and nCa*100 == 1
print('INO-ALK-019 OK')
