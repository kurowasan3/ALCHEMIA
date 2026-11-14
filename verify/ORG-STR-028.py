# 検証スクリプト ORG-STR-028  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-STR-028/item.json'
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
# ---- 個別の検証 ORG-STR-028 ----
import math
R = 8.31; muG = 2.60
def xg(dE, T):
    r = math.exp(-dE/(R*T)); return 2*r/(1+2*r)
# (1)
q = 4990/(R*300); print('ΔE/RT', q); assert abs(q-2.00) < 0.005
x1 = 2*0.135/(1+2*0.135); print('x1', x1, xg(4990,300)); assert sig(x1,3) == 0.213 and sig(xg(4990,300),3) == 0.213
mu1 = muG*0.462; assert sig(mu1,3) == 1.20; assert sig(muG*math.sqrt(x1),3) == 1.20
assert abs(math.sqrt(0.213)-0.462) < 0.0006
# (2)
assert abs(4990/(R*600) - 1.00) < 0.003
x2 = 2*0.368/(1+2*0.368); print('x2', x2, xg(4990,600)); assert sig(x2,3) == 0.424 and sig(xg(4990,600),3) == 0.424
assert abs(math.sqrt(0.424)-0.651) < 0.0006
mu2 = muG*0.651; print('mu2', mu2); assert sig(mu2,3) == 1.69 and mu2 > mu1
# (3)
assert abs(xg(4990, 1e9) - 2/3) < 1e-4
# (4)
x4 = (1.84/2.60)**2; print('x4', x4); assert sig(x4,2) == 0.50
dE4 = R*300*0.693; print('dE4', dE4); assert sig(dE4/1000, 2) == 1.7
assert abs(xg(R*300*math.log(2), 300) - 0.5) < 1e-12
print('ALL OK')
