# 検証スクリプト THE-MOL-029  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-MOL-029/item.json'
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
x, y = sp.symbols('x y')
assert sp.solve(sp.Eq(sp.Rational(20,100)*x + 1*(sp.Rational(4,10)-x), sp.Rational(1,2)*sp.Rational(4,10)), x) == [sp.Rational(1,4)]
assert sp.solve(sp.Eq(sp.Rational(20,100)*x, sp.Rational(1,10)*sp.Rational(4,10)), x) == [sp.Rational(1,5)]
assert sp.solve(sp.Eq(1*x, sp.Rational(1,10)*sp.Rational(4,10)), x) == [sp.Rational(1,25)]
assert sp.solve(sp.Eq(sp.Rational(1,10)*y + sp.Rational(3,10)*(200-y), sp.Rational(1,4)*200), y) == [50]
assert sp.Rational(250,150) == sp.Rational(5,3) and sp.Rational(50,150) == sp.Rational(1,3)

print('ALL OK')
