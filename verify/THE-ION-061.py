# 検証スクリプト THE-ION-061  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ION-061/item.json'
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
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(re.findall(r'\((\d+)\)', it['answer'])) == len(nums), '答え中に余分な (n)'
    if n_items: assert nums == [str(i+1) for i in range(n_items)], (n_items, nums)
    else: assert not nums
    print('小問数', n_items, '答え', nums)
check_structure()
check_all_equations()
def bis(f, lo=-16, hi=1):
    """f(h) が h について増加 → f=0 の h（log 空間の二分法）"""
    for _ in range(200):
        m = (lo + hi) / 2
        if f(10**m) > 0: hi = m
        else: lo = m
    return 10**((lo + hi) / 2)
# ---- 個別の検証 ----
from math import log10
Kb=2.0e-5; Kw=1.0e-14; Ka=Kw/Kb
def sol(CN, CCl):
    r = bis(lambda h: CN*h/(h+Ka) + h - CCl - Kw/h)
    nh4 = CN*r/(r+Ka)
    return {'NH4+':nh4, 'Cl-':CCl, 'NH3':CN-nh4, 'H+':r, 'OH-':Kw/r}
R = sol(2.0/30, 1.0/30); P = sol(2.0/40, 2.0/40); Q = sol(2.0/50, 3.0/50)
order = lambda d: sorted(d, key=lambda k: -d[k])
print('R', order(R)); assert order(R) == ['NH4+','Cl-','NH3','OH-','H+']
assert abs(-log10(R['H+']) - 9.30) < 0.01
assert round(-log10(P['H+']), 2) == 5.30
assert P['H+'] > P['NH3'] and abs(P['H+'] - P['NH3'] - P['OH-']) < 1e-12
assert round(-log10(Q['H+']), 2) == 1.70
assert sig(Q['NH3'], 2) == 1.0e-9
assert order(Q) == ['Cl-','NH4+','H+','NH3','OH-']
# 近似値（解説）
assert abs((5.0e-10*0.050)**0.5 - 5.0e-6) < 1e-12

print('ALL OK')
