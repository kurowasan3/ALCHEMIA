# 検証スクリプト ORG-HYD-028  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-HYD-028/item.json'
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
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ORG-HYD-028 ----
import sympy as sp
y, z = sp.symbols('y z')
x = sp.Rational(240,100)/240   # Ag2C2 2.40 g / 240
assert 108*2 + 12*2 == 240
sol = sp.solve([x + y + z - sp.Rational(1,10), y + 2*z - sp.Rational(15,100)], [y, z])
print('C2H2', x, 'C2H4', sol[y], 'C2H6', sol[z])
assert x == sp.Rational(1,100) and sol[y] == sp.Rational(3,100) and sol[z] == sp.Rational(6,100)
br = 2*x + sol[y]; assert br == sp.Rational(5,100)
assert sp.Rational(1,10)/sp.Rational(25,100) == sp.Rational(2,5)
# (5) H2:C2H2=1:1 のとき 2x+y が一定であることを記号的に確認
a, b, c = sp.symbols('a b c', nonnegative=True)
s2 = sp.solve([a+b+c-1, b+2*c-1], [a, b])
assert sp.simplify(2*s2[a] + s2[b]) == 1
print('ALL OK')
