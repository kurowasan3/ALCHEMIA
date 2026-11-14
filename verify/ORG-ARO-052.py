# 検証スクリプト ORG-ARO-052  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-ARO-052/item.json'
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
from rdkit import Chem
from itertools import product
from fractions import Fraction as F
import sympy as sp
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
def ring(pos):
    parts = []
    for i in range(6):
        a = 'c' + ('1' if i in (0, 5) else '')
        if i in pos: a += '(C)'
        parts.append(a)
    return can(''.join(parts))
NAMES = {can('c1ccccc1'):'benzene', can('Cc1ccccc1'):'toluene', can('Cc1ccccc1C'):'o', can('Cc1cccc(C)c1'):'m',
         can('Cc1ccc(C)cc1'):'p', can('Cc1ccc(C)c(C)c1'):'124', can('Cc1cc(C)cc(C)c1'):'135'}
def dist(x, r):
    out = {}
    for kinds in product('AP', repeat=3):
        pk = F(1)
        for k in kinds: pk *= x if k == 'P' else 1 - x
        opts = [[None] if k == 'A' else [2*u, 2*u+1] for u, k in enumerate(kinds)]
        ws = {}; tot = F(0)
        for ch in product(*opts):
            pos = {c for c in ch if c is not None}
            mm = sum(1 for a, b in ((1, 2), (3, 4), (5, 0)) if a in pos and b in pos)
            w = r**mm; tot += w
            n = NAMES[ring(pos)]; ws[n] = ws.get(n, 0) + w
        for n, w in ws.items(): out[n] = out.get(n, 0) + pk*w/tot
    return out
d = dist(F(1), F(1)); assert d['124'] == F(3,4) and d['135'] == F(1,4)
for x in (F(1,5), F(1,2), F(2,3)):
    d = dist(x, F(1))
    exp = {'benzene':(1-x)**3, 'toluene':3*x*(1-x)**2, 'o':F(3,4)*x**2*(1-x), 'm':F(3,2)*x**2*(1-x),
           'p':F(3,4)*x**2*(1-x), '124':F(3,4)*x**3, '135':F(1,4)*x**3}
    assert all(d.get(k, 0) == v for k, v in exp.items()) and sum(d.values()) == 1
X = sp.symbols('x'); assert set(sp.solve(sp.Rational(3,2)*X**2*(1-X) - 3*X*(1-X)**2, X)) == {0, sp.Rational(2,3), 1}
d = dist(F(2,3), F(1)); assert d['m'] == d['toluene'] == F(2,9); assert d['o'] + d['124'] == F(1,3)
d = dist(F(1), F(1,2)); assert d['124']/d['135'] == F(3,2)
d = dist(F(1,2), F(1,2)); assert (d['o'], d['m'], d['p']) == tuple(F(k,7)*(d['o']+d['m']+d['p']) for k in (1,4,2))
d = dist(F(1,3), F(1,2)); assert d['o']*4 == d['m'] and d['p']*2 == d['m']
for e in ['3C2H2 -> C6H6', '3C3H4 -> C9H12']: assert check_equation(e) is None
print('生成比のモデル計算 OK')

print('ALL OK')
