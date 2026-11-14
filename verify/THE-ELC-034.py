# 検証スクリプト THE-ELC-034  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-ELC-034/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    return re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', s)
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq
def should_check(eq):
    # 矢印が1つの完全な式はすべて（validate ツールが一般式とみなして飛ばす Zn・Sn なども含めて）検査する．
    # 一般式・同位体の式は矢印に空ラベル [\\text{}] を付けてあるので除く．
    if re.search(r'boxed|□|\?|？|…|\\square', eq): return False
    if re.search(r'(->|<=>|<-)\[\\text\{\}\]', eq): return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1: return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None
def check_all_equations():
    eqs = sorted({b for k in ('problem','answer','solution') for b in ce_blocks(it[k]) if should_check(b)})
    for e in eqs:
        r = check_equation(clean_ce(e)); print('式', e, '→', r or 'OK'); assert r is None, e
    return eqs
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items+1)), nums
    print('小問数', n_items)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
EQS = check_all_equations()

# ---- 個別の検証 ----

from sympy import Rational as R, symbols, solve
assert '2H2O + 2e- -> H2 + 2OH-' in EQS and '2Cl- -> Cl2 + 2e-' in EQS
w = symbols('w')
ws = solve(40 - R(471,1000)*(40 + 18*(w-1)), w)[0]; print('w', float(ws)); assert sig(float(ws),3) == 3.50
w = R(35,10); n = 2
f = 80/R(32,100) - 80 - 18*(w-1)*n; print('feed', f); assert f == 80
f0 = 80/R(32,100) - 80 + 18*n; assert f0 == 206
real = R(80, 250 + 18*(w-1)*n + 0) ; # 給水 206 g のときの実際の濃度
real = R(80, 80 + 206 + 18*(w-1)*n); print('実際', float(real)); assert sig(float(real)*100,3) == 21.3
nacl, h2o = 58.5*n, 18*w*n; assert nacl == 117 and h2o == 126
c3 = R(260-117, (260-117)+(740-126)); print('c3', float(c3)); assert sig(float(c3)*100,3) == 18.9
assert sig(143/883*100,3) == 16.2
x = symbols('x'); xs = solve(R(26,100)*x - 117 - R(2,10)*(x - 243), x)[0]; print('x', xs); assert xs == 1140
wt = symbols('wt'); wts = solve(R(585,10) - R(26,100)*(R(585,10) + 18*wt), wt)[0]; print('w閾値', float(wts)); assert sig(float(wts),3) == 9.25
# 上限の確認：f>0 なら濃度は 47.1% 未満
for fv in (0, 10, 100): print(fv, float(R(80, 80 + 18*(R(3496,1000)-1)*2 + fv)))

print('全検証OK:', ITEM)
