ITEM_ID = 'THE-MOL-022'
# ---- 共通ヘッダ（C0xx 検証スクリプト） ----
import sys, json, re, os
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational, nsimplify, N, symbols, solve, Eq
ITEM = os.path.join('/home/claude/chem/out/items', ITEM_ID, 'item.json')
it = json.load(open(ITEM, encoding='utf-8'))
FAIL = []
def ok(cond, msg):
    print(('OK   ' if cond else 'FAIL ') + msg)
    if not cond: FAIL.append(msg)
def near(name, val, exp, rel=0.006):
    ok(abs(float(val) - float(exp)) <= rel * abs(float(exp)), f'{name}: computed {float(val):.6g} vs stated {exp}')
def sig(x, n):
    from math import log10, floor
    x = float(x); return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def ce_blocks(s):
    out = []; i = 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k = j + 4; d = 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j + 4:k - 1]); i = k
def check_equations():
    txt = ' '.join(it[k] for k in ('problem', 'answer', 'solution'))
    n = 0
    for c in ce_blocks(txt):
        if c.count('->') + c.count('<=>') != 1 or '\\boxed' in c or re.search(r'\bR\b|\bn\b|\bM\b|\bx\b', c):
            continue
        r = check_equation(c.replace('<=>', '->')); n += 1
        ok(r is None, f'equation [{c}] ' + ('' if r is None else str(r)))
    print(f'  ({n} equations checked)')
def check_parts():
    p = it['problem']; a = it['answer']
    ni = p.count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', a)
    ok(p.count('\\begin{enumerate}') <= 1, 'enumerate not nested / single')
    if ni:
        ok([int(x) for x in nums] == list(range(1, ni + 1)), f'answer numbering {nums} matches {ni} items')
    else:
        ok(not nums, 'single answer without numbering')
    for k in [] or ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        if k not in it: ok(False, f'key {k}')
    ok(it['difficulty'] in (1, 2, 3), 'difficulty 1-3')
    ok(all(f in ('qa','equation','compare','condition','pathway','calc','structure') for f in it['format'].split(',')) and len(it['format'].split(','))<=2, 'format valid')
    for f in ('problem_fig', 'solution_fig'):
        if it[f]: ok(os.path.exists(os.path.join(os.path.dirname(ITEM), it[f])), f'{f} exists')
def finish():
    check_equations(); check_parts()
    print('RESULT:', 'ALL OK' if not FAIL else f'{len(FAIL)} FAIL')
    sys.exit(1 if FAIL else 0)
# ---- 個別の検証 ----
for name,(a,b) in {'A':(2,1),'B':(1,1),'C':(1,2)}.items():
    wN = Rational(14*a, 14*a+16*b)*100
    print(name, 'N%', float(wN), 'O per gN', float(Rational(16*b,14*a)))
near('A N%', Rational(28,44)*100, 63.64, rel=3e-4); near('B N%', Rational(14,30)*100, 46.67, rel=3e-4); near('C N%', Rational(14,46)*100, 30.43, rel=3e-4)
near('A O/gN', Rational(16,28), 0.571, rel=0.002); near('B', Rational(16,14), 1.14, rel=0.004); near('C', Rational(32,14), 2.29, rel=0.002)
near('A from %', Rational(3636,6364), 0.571, rel=0.001)
near('B from %', Rational(5333,4667), 1.14, rel=0.004); near('C from %', Rational(6957,3043), 2.29, rel=0.002)
ok(Rational(16,28)*2 == Rational(16,14) and Rational(16,14)*2 == Rational(32,14), 'ratio 1:2:4')

finish()
