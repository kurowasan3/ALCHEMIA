# 検証スクリプト POL-PRT-061  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/POL-PRT-061/item.json'
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
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


from fractions import Fraction
PAIR = {'A':'T','T':'A','G':'C','C':'G','U':'A','H':'C'}
def replicate(mols):
    out = []
    for a, b in mols:
        out.append((a, PAIR[a])); out.append((b, PAIR[b]))
    return out
def run(pair, n):
    mols = [pair]
    for _ in range(n): mols = replicate(mols)
    return mols
norm = lambda p: tuple(sorted(p))
m2 = run(('G','U'), 2); print(m2)
assert sum(norm(p)==('A','T') for p in m2) == 1
m3 = run(('G','U'), 3); f = Fraction(sum(norm(p)==('A','T') for p in m3), len(m3)); print(f); assert f == Fraction(3,8)
for n in range(1,7):
    assert sum(norm(p)==('A','T') for p in run(('G','U'), n)) == 2**(n-1)-1
m4 = run(('H','T'), 2); print(m4); assert sum(norm(p)==('C','G') for p in m4) == 1
print('OK')
