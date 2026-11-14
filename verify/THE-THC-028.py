# 検証スクリプト THE-THC-028  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-THC-028/item.json'
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
    keys = ['src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig']
    assert all(k in it for k in keys), [k for k in keys if k not in it]
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

from sympy import Rational as R, symbols, solve
c = R(42,10); C = symbols('C')
Cc = solve(100*c*(R(400,10)-R(296,10)) - 100*c*(R(296,10)-20) - C*(R(296,10)-20), C)[0]
assert Cc == 35
Q = (200*c + Cc)*R(66,10); assert Q == 5775
dH = Q/R(1,10)/1000; assert sig(float(dH),2) == 58
dH0 = 200*c*R(66,10)/R(1,10)/1000; assert sig(float(dH0),2) == 55
dT = dH*1000*R(5,100)/(150*c + Cc); print(float(dT)); assert sig(float(dT),2) == 4.3
assert check_equation('HCl + NaOH -> NaCl + H2O') is None
print('J144 OK')
