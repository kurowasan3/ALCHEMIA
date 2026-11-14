ITEM_ID = 'THE-GAS-021'
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
R=8300
ok(sig(Rational(4,100)*R*300/2,2)==50000,'(1)')
ok(sig(Rational(4,100)*R*300/5,2)==20000,'(2)')
P = Rational(24,100)*R/(Rational(2,300)+Rational(3,500)); nA = P*2/(R*300)
ok(sig(P,2)==1.6e5, f'(3) P {float(P):.0f}'); ok(sig(nA,2)==0.13, f'(3) nA {float(nA):.4f}'); ok(nA/Rational(24,100) > Rational(2,5), '(3) larger')
P4 = Rational(24,100)*R*500/5; ok(sig(P4,2)==2.0e5, f'(4) {float(P4):.0f}')
ng = Rational(3600*3)/(R*300); liq = Rational(8,100)*Rational(3,5) - ng; ok(sig(liq,2)==0.044, f'(5) {float(liq):.5f}')
ok(Rational(20,100) >= 2*Rational(4,100), 'O2 excess')

finish()
