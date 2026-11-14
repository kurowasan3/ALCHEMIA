# 検証スクリプト THE-GAS-031  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-GAS-031/item.json'
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
assert check_equation('CH3OCH3 -> CH4 + CO + H2') is None
n = 4.0e4*2.0/(8.31e3*400); print('n', n); assert sig(n, 2) == 2.4e-2
P0 = Q(4)*10**4 * Q(800, 400); assert P0 == 8*10**4
a = symbols('a'); al = solve(P0*(1 + 2*a) - Q(144)*10**3, a)[0]; print('alpha', al); assert al == Q(2, 5)
pCH4 = P0*al; assert pCH4 == 32000
M = {'C2H6O': 2*12 + 6*1.0 + 16, 'CH4': 12 + 4.0, 'CO': 28, 'H2': 2.0}
Mavg = ((1-float(al))*M['C2H6O'] + float(al)*(M['CH4'] + M['CO'] + M['H2'])) / (1 + 2*float(al))
print('平均分子量', Mavg); assert abs(Mavg - 46/1.8) < 1e-9 and sig(Mavg, 2) == 26
# 誤り例：温度補正なしだと alpha = (1.44e5/4.0e4 - 1)/2 = 1.3 (>1，不可能)
print('補正なしの alpha', (1.44e5/4.0e4 - 1)/2)

check_equations()
print('全検証OK:', ITEM)
