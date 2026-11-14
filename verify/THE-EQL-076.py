# 検証スクリプト THE-EQL-076  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/THE-EQL-076/item.json'
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
S, S0, V, K, t = sp.symbols('S S0 V K t', positive=True)
# 式②を時間で微分すると -dS/dt が式①になる
tS = (S0 - S)/V + K/V*sp.log(S0/S)
assert sp.simplify(-1/sp.diff(tS, S) - V*S/(K+S)) == 0
f = lambda s0, s, v, k: (s0 - s)/v + k/v*sp.log(sp.Rational(s0)/sp.Rational(s))
v, k = sp.Rational(1,10), 2
approx = (20 - 2)/v + sp.log(10)/(v/k)
assert sig(float(approx), 2) == 230, float(approx)
exact = f(20, sp.Rational(1,5), v, k)
assert sig(float(exact), 2) == 290, float(exact)
seg1, seg2 = f(20, 2, v, k), f(2, sp.Rational(1,5), v, k)
assert float(seg1) > 180 and float(seg2) > float(sp.log(10)/(v/k))
print('近似', float(approx), '正確', float(exact), '区間', float(seg1), float(seg2))
t2E = f(20, sp.Rational(1,5), 2*v, k); t2S = f(40, sp.Rational(1,5), v, k)
assert sig(float(t2E), 2) == 150 and sig(float(t2S), 2) == 500, (float(t2E), float(t2S))
# 与えた対数値でも同じ丸め
ln2, ln10 = 0.693, 2.303
assert sig(198 + 20*2*ln10, 2) == 290 and sig(398 + 20*(2*ln10+ln2), 2) == 500
# (5)
th = lambda s0, v, k: s0/(2*v) + k*ln2/v
print('t1/2', th(2.0, 0.25, 2.0), th(10, 0.25, 2.0))
assert round(th(2.0,0.25,2.0),1) == 9.5 and round(th(10,0.25,2.0),1) == 25.5
slope = (25.5-9.5)/(10-2.0); V5 = 1/(2*slope); icpt = 9.5 - slope*2.0; K5 = icpt*V5/ln2
assert abs(V5-0.25) < 1e-12 and sig(K5, 2) == 2.0, (V5, K5)
# 真の速さは近似より小さい
for s in [0.1, 1, 2, 10, 20]:
    assert 0.1*s/(2+s) < 0.1 and 0.1*s/(2+s) < 0.1*s/2
print('ALL checks for enzyme kinetics OK')

print('ALL OK')
