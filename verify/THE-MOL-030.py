# 検証スクリプト THE-MOL-030  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-MOL-030/item.json'
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
nH2max = F(224,1000)/F(224,10)
assert nH2max == F(1,100)
assert F(24,100)/24 == nH2max
c = nH2max / F(50,1000); assert c == F(1,5)
# (2) 25 mL
assert nH2max/2*22400 == 112 and nH2max/2*24 == F(12,100)
# (3) HCl 0.20 mol/L 50 mL -> H+ 0.010 mol -> H2 0.005
nHCl = F(1,5)*F(50,1000); assert nHCl/2*22400 == 112 and nHCl/2*24 == F(12,100)
# (4) Al
assert (F(1,27)*F(3,2))/(F(1,24)) == F(4,3)
assert nH2max*F(2,3)*27 == F(18,100)
# graph points
pts = [min(F(6,100)*k/24, nH2max)*22400 for k in range(1,8)]
assert pts == [56,112,168,224,224,224,224]
for e in ['Mg + 2HCl -> MgCl2 + H2']:
    assert check_equation(e) is None
import os
assert os.path.exists(os.path.join(os.path.dirname(ITEM), 'problem.svg'))

print('ALL OK')
