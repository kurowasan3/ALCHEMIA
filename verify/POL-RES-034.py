# 検証スクリプト POL-RES-034  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/POL-RES-034/item.json'
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|\$', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式の検査', n, '本 OK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

import itertools
from collections import Counter
from sympy import symbols, solve, Rational
def canon(s): return min(s, s[::-1])
# (1)
tri = {canon(''.join(t)) for t in itertools.product('mr', repeat=2)}
pen = {canon(''.join(t)) for t in itertools.product('mr', repeat=4)}
print('(1)', sorted(tri), len(pen)); assert len(tri) == 3 and len(pen) == 10
# (2)
p = symbols('p', positive=True)
pv = [s for s in solve(p**2 - Rational(81, 100), p) if s > 0][0]
assert pv == Rational(9, 10) and (1 - pv)**2 == Rational(1, 100)
assert sig(float(pv**4), 2) == 0.66
# (3)
def dyads(sides): return ''.join('m' if sides[i] == sides[i+1] else 'r' for i in range(len(sides)-1))
def windows(seq): return Counter(canon(seq[i:i+4]) for i in range(len(seq)-3))
a = dyads('FFFFFBFFFFF'); print('(3a)', a, windows(a))
wa = windows(a); assert wa['mmmr'] == 2 and wa['mmrr'] == 2 and wa['mrrm'] == 1 and 'mmrm' not in wa
b = dyads('FFFFFBBBBBB'); print('(3b)', b, windows(b))
wb = windows(b); assert wb['mmmr'] == 2 and wb['mmrm'] == 2 and 'mmrr' not in wb
assert '\\cdots mmrrmm\\cdots' in it['answer'] and '\\cdots mmrmm\\cdots' in it['answer']
# (4) 誤り率 2% の孤立欠陥（50 単位ごとに 1 個）
sides = ''.join('B' if i % 50 == 25 else 'F' for i in range(50*400))
dy = dyads(sides)
# 単位 j を中心とするペンタッド：単位 j-2..j+2 → 記号 j-2..j+1
cnt = Counter(canon(dy[j-2:j+2]) for j in range(2, len(sides)-2))
tot = sum(cnt.values())
fr = {k: v/tot for k, v in cnt.items()}
print('(4)', {k: round(v, 4) for k, v in fr.items()})
assert abs(fr['mmmm'] - 0.90) < 1e-3 and abs(fr['mmmr'] - 0.04) < 1e-3 and abs(fr['mmrr'] - 0.04) < 1e-3 and abs(fr['mrrm'] - 0.02) < 1e-3
# (5)
P = 0.95
prob = Counter()
for t in itertools.product('mr', repeat=4):
    pr = 1
    for c in t: pr *= P if c == 'm' else 1 - P
    prob[canon(''.join(t))] += pr
print('(5)', prob['mmrm'], prob['mmrr'], prob['mmrm']/prob['mmrr'])
assert sig(prob['mmrm'], 2) == 0.086 and sig(prob['mmrr'], 2) == 0.0045

print('ALL OK')
