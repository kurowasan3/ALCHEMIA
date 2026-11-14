# 検証スクリプト POL-RES-018  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/POL-RES-018/item.json'
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
from fractions import Fraction as F
import sympy as sp
w = F(1); M1, M2 = F(10**4), F(10**5)
N1, N2 = w/M1, w/M2
Mn = (N1*M1 + N2*M2)/(N1 + N2); Mw = (N1*M1**2 + N2*M2**2)/(N1*M1 + N2*M2)
assert round(float(Mn), -2) == 18200 and Mw == 55000 and Mw/Mn == F(121,40)
ws = [(F(495,10), M1), (F(495,10), M2), (F(1), F(1000))]
Mn2 = sum(a for a, m in ws)/sum(a/m for a, m in ws); Mw2 = sum(a*m for a, m in ws)/sum(a for a, m in ws)
assert round(float(Mn2), -2) == 15500 and Mw2 == 54460
assert round(float((F(1)/1000)/sum(a/m for a, m in ws)), 3) == 0.155
p, n = sp.symbols('p n', positive=True)
ser1 = sp.series(1/(1-p)**2, p, 0, 30).removeO(); ser2 = sp.series((1+p)/(1-p)**3, p, 0, 30).removeO()
for k in range(1, 30):
    assert ser1.coeff(p, k-1) == k and ser2.coeff(p, k-1) == k**2
M0 = 113; pp = 0.99
Mn_num = sum(n_*M0*(1-pp)*pp**(n_-1) for n_ in range(1, 20000))
Mw_num = sum((n_*M0)**2*(1-pp)*pp**(n_-1) for n_ in range(1, 20000))/Mn_num
assert abs(Mn_num - 11300) < 1e-6 and abs(Mw_num - 113*1.99/0.01) < 1e-6
assert round(113*1.99/0.01, -2) == 22500
P = F(99, 100); f = lambda k: k*P**(k-1)
assert f(99) == f(100) and f(98) < f(99) and f(101) < f(100)
print('平均分子量 OK')
assert round((1 - Mn2/Mn)*100) == 15 and round((1 - Mw2/Mw)*100) == 1

print('ALL OK')
