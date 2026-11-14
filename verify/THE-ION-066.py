# 検証スクリプト THE-ION-066  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/THE-ION-066/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, formula_mass
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
from sympy import symbols, nsolve, Rational, sqrt, solve, Float
Ka = Rational(7, 10000); K2 = 4
h, a = symbols('h a', positive=True)
def exact(C, g):
    f = Ka*a/h; b = K2*a*f
    s = nsolve([f + b - h, a + f + 2*b - C], [h, a], g)
    return float(s[0]), float(s[1])
# (1)
x = symbols('x', positive=True)
h1 = float([s for s in solve(x**2 + Ka*x - Ka, x) if s > 0][0]); print('(1)', h1)
assert sig(h1, 2) == 0.026 and sig(float(sqrt(Ka)), 2) == 0.026
# (2) 厳密解が 0.0540 になること，逆算
hE, aE = exact(1, [0.05, 0.9]); print('exact h', hE); assert sig(hE, 3) == 0.0540
H = 0.0540; A = (1 - 2*H)/(1 - float(Ka)/H); Fm = float(Ka)*A/H; B = H - Fm; k2 = B/(A*Fm)
print('(2)', A, Fm, B, k2)
assert (sig(A, 2), sig(Fm, 2), sig(B, 2), sig(k2, 2)) == (0.90, 0.012, 0.042, 4.0)
assert sig(B/H, 2) == 0.78
# (3) 境界
fb = (float(Ka)*0.0125/1.05)**0.5; Cb = 0.0125 + fb + 2*0.05*fb; print('(3)', fb, Cb)
hB, aB = exact(Cb, [0.003, 0.012]); fB = float(Ka)*aB/hB; print('  check ratio', K2*aB)
assert abs(K2*aB - 0.05) < 1e-3 and sig(Cb, 2) == 0.016  # a=C と誤ると 0.0125（四捨五入で 1.3e-2）
# (4)
h4, a4 = exact(Rational(1, 10), [0.01, 0.09]); print('(4)', h4, a4, hE/h4)
assert sig(h4, 2) == 0.0092 and sig(hE/h4, 2) == 5.9
h4n = float([s for s in solve(x**2 + Ka*x - Ka*Rational(1, 10), x) if s > 0][0])
print('  no HF2-', h4n, h1/h4n); assert sig(h1/h4n, 2) == 3.3
# 反復
aa = 0.100
for i in range(3):
    hh = (float(Ka)*aa*(1 + K2*aa))**0.5; ff = float(Ka)*aa/hh; bb = K2*aa*ff; print('  iter', aa, hh, ff, bb); aa = 0.100 - ff - 2*bb
assert abs(1 + K2*aE - 4.6) < 0.05 and abs(1 + K2*a4 - 1.35) < 0.01
for s in ['2.6\\times10^{-2}', '0.90', '1.2\\times10^{-2}', '4.2\\times10^{-2}', '4.0', '1.6\\times10^{-2}', '9.2\\times10^{-3}', '1/5.9', '1/3.3']:
    assert s in it['answer'], s

print('ALL OK')
