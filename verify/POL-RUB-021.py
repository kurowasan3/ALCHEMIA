# 検証スクリプト POL-RUB-021  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/POL-RUB-021/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
import sympy as sp
x, y, h, V = sp.symbols('x y h V', positive=True)
K = sp.Rational(3,2)
expr = (x*(x/V))/((2-x)*((2-x)/V))
assert sp.simplify(sp.diff(expr, V)) == 0
x1 = [s for s in sp.solve(sp.Eq(x**2, K*(2-x)**2), x) if 0 < s < 2][0]
assert sig(float(x1/2*100), 2) == 55
assert abs(2*1.22/2.22 - 1.10) < 0.005
x1v = float(x1)
y1 = [s for s in sp.solve(sp.Eq(y*(x1v+y), 1.5*(2-y)*(2-x1v-y)), y) if 0 < s < 2-x1v][0]
tot = (x1v + float(y1))/2*100
assert sig(tot, 2) == 81, tot
assert abs((5.45-4.93) - 0.52) < 1e-9 and abs(0.5*0.52**2 - 5.45*0.52 + 2.70) < 0.01
hs = sp.solve(sp.Eq(sp.Rational(2,100)*(h-sp.Rational(198,100)), K*sp.Rational(198,100)**2), h)[0]
assert sig(float(hs), 2) == 300, hs
print('イオン交換平衡 OK', float(x1), float(y1), tot, float(hs))

print('ALL OK')
