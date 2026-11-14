# 検証スクリプト THE-GAS-030  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-GAS-030/item.json'
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡'.split()
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
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and all(x.strip() for x in re.split(r'->|<=>|<-', b)):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    assert it['field'] in FIELDS
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    for f in ('problem_fig', 'solution_fig'):
        if it[f]:
            import os
            p = os.path.join(os.path.dirname(ITEM), it[f]); assert os.path.exists(p), p
            t = open(p).read(); assert 'viewBox' in t and '<script' not in t and 'foreignObject' not in t
    print('小問数', n_items, '答え数', n_ans)
check_structure()
from sympy import symbols, solve, Rational as Q
Pex = Q(10)**5; VI = Q(2)
def final(PI, PII, VII):
    PV = PI*VI + PII*VII
    Vtot = PV/Pex
    if Vtot >= VI: return Pex, Vtot - VI
    return PV/VI, 0
P, V = final(2*Pex, Pex, 3); print('(1)', P, V); assert P == Pex and V == 5 and V > 3
nN, nA = 2*Pex*VI, Pex*3
pN = Pex*nN/(nN+nA); pA = Pex*nA/(nN+nA); print('(2)', float(pN), float(pA))
assert sig(float(pN), 2) == 5.7e4 and sig(float(pA), 2) == 4.3e4
p = symbols('p'); ps = solve(p*VI + Pex*3 - Pex*(VI + Q(12, 10)), p)[0]; print('(3)', ps)
assert ps == 10**4 and final(ps, Pex, 3) == (Pex, Q(12, 10))
P4, V4 = final(0, Pex, Q(3, 2)); print('(4)', float(P4), V4); assert V4 == 0 and sig(float(P4), 2) == 7.5e4

check_equations()
print('全検証OK:', ITEM)
