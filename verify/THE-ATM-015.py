# 検証スクリプト THE-ATM-015  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-ATM-015/item.json'
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
def nuc(tok):
    tok = tok.strip()
    if tok == 'e-': return (0, -1)
    m = re.fullmatch(r'\^\{(\d+)\}_\{(\d+)\}([A-Za-z]+)', tok); assert m, tok
    return (int(m.group(1)), int(m.group(2)))
nuc_eqs = [b for k in ('problem','answer','solution') for b in ce_blocks(it[k]) if '->' in b]
for e in nuc_eqs:
    lhs, rhs = re.split(r'->\[\\text\{\}\]', e)
    L = [nuc(t) for t in lhs.split(' + ')]; R = [nuc(t) for t in rhs.split(' + ')]
    assert sum(a for a,_ in L) == sum(a for a,_ in R) and sum(z for _,z in L) == sum(z for _,z in R), e
assert len(nuc_eqs) == 2; print('核反応式 OK', nuc_eqs)
# (2) の式を時間発展の数値積分（微分方程式）で独立に確かめる
t, T, N0 = sp.symbols('t T N0', positive=True)
lam = sp.log(2)/T
N = N0*sp.exp(-lam*t)
Ar = sp.integrate(sp.Rational(107,1000)*lam*N, (t, 0, t))
Rexpr = sp.simplify(Ar/N)
assert sp.simplify(Rexpr - sp.Rational(107,1000)*(2**(t/T)-1)) == 0
Tv = 1.25e9
def age(R, f=0.107): return Tv*sp.log(1+R/f, 2)
tA = float(age(0.321)); print('A', tA); assert sig(tA,2) == 2.5e9
tB = float(age(0.321*0.70)); print('B exact', tB)
tB_tab = Tv*0.491/0.301; print('B tab', tB_tab); assert sig(tB,2) == sig(tB_tab,2) == 2.0e9
tC = float(age(0.321, 1.0)); tC_tab = Tv*0.121/0.301; print('C', tC, tC_tab); assert sig(tC,2) == sig(tC_tab,2) == 5.0e8
# 与えた対数値の確認
import math
assert abs(math.log10(2)-0.301) < 5e-4 and abs(math.log10(3.1)-0.491) < 5e-4 and abs(math.log10(1.32)-0.121) < 5e-4

print('ALL OK')
