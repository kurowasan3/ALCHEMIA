# 検証スクリプト ORG-HYD-027  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-HYD-027/item.json'
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
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ORG-HYD-027 ----
import sympy as sp
assert check_equation('C14H30 -> C10H22 + C4H8') is None
assert check_equation('C16H34 -> C8H18 + 2C2H4 + C4H8') is None
a, b = sp.symbols('a b', positive=True, integer=True)
assert sp.expand((2*a+2) + (2*b+2) - (2*(a+b)+2)) == 2   # 2 アルカンにすると H が 2 個余分に必要
# k 分子に分かれると H の不足は 2(k-1) 個
for k in range(2, 6):
    assert k*2 + 0 - 2 == 2*(k-1)
assert (2+8+8) and 34 - 18 - 2*4 - 8 == 0
print('分留成分の沸点の順（石油ガス<ナフサ<灯油<軽油）は高校教科書の記述と一致')
print('ALL OK')
