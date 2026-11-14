ITEM_ID = 'THE-GAS-017'
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
PB = (Rational(25,10)*10**5*8 - 10**5*3)/5; near('(1)', PB, 3.4e5)
Ar = (Rational(1,2)*10**5*3 + Rational(1,4)*PB*5)/8; near('(2)', Ar, 7.19e4, rel=0.001); ok(sig(Ar,2)==7.2e4,'(2) sig')
RT=8300*300
nA = Rational(10**5*3, RT); nB = PB*5/RT; nCH4 = nA/2; nO2 = Rational(3,4)*nB
ok(nO2 > 2*nCH4, 'O2 in excess')
ntot = nA + nB  # unchanged by combustion
P3 = ntot*8300*500/8; near('(3)', P3, 4.17e5, rel=0.002); ok(sig(P3,2)==4.2e5,'(3) sig')
P4 = (nO2-2*nCH4)*8300*500/8; near('(4)', P4, 2.03e5, rel=0.002); ok(sig(P4,2)==2.0e5,'(4) sig')

finish()
