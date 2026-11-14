# 検証スクリプト THE-ELC-036  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ELC-036/item.json'
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
Q = 0.500*(32*60+10); assert Q == 965.0
ne = 2*0.318/63.5; F = Q/ne; print(F); assert sig(F,3) == 9.63e4
NA = sig(F,3)/1.60e-19; print(NA); assert sig(NA,3) == 6.02e23
assert sig(F/1.60e-19,3) == 6.02e23
# (a) 水素が出ると析出量が減る→F 大
assert Q/(2*0.300/63.5) > F
# (b) ぬれていると質量増加が大→F 小
assert Q/(2*0.330/63.5) < F

print('ALL OK')
