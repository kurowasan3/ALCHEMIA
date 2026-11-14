# 検証スクリプト ORG-ALC-013  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ALC-013/item.json'
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
# ---- 個別の検証 ORG-ALC-013 ----
from collections import Counter
for e in ['CH2=CH2 + PdCl2 + H2O -> CH3CHO + Pd + 2HCl', 'Pd + 2CuCl2 -> PdCl2 + 2CuCl', '4CuCl + 4HCl + O2 -> 4CuCl2 + 2H2O',
          '2CH2=CH2 + O2 -> 2CH3CHO', '2CH2=CHCH3 + O2 -> 2CH3COCH3']:
    assert check_equation(e) is None, e
# 段階 1×2 + 段階 2×2 + 段階 3 で全体の式になることを，化学種の数の足し算で確認
def side(s): 
    c = Counter()
    for t in s.split(' + '):
        import re
        m = re.match(r'(\d*)(.+)', t.strip()); c[m.group(2)] += int(m.group(1) or 1)
    return c
steps = [('CH2=CH2 + PdCl2 + H2O', 'CH3CHO + Pd + 2HCl', 2), ('Pd + 2CuCl2', 'PdCl2 + 2CuCl', 2), ('4CuCl + 4HCl + O2', '4CuCl2 + 2H2O', 1)]
L, R = Counter(), Counter()
for l, r, k in steps:
    for sp_, n in side(l).items(): L[sp_] += k*n
    for sp_, n in side(r).items(): R[sp_] += k*n
net = {s: R[s]-L[s] for s in set(L)|set(R) if R[s] != L[s]}
print(net); assert net == {'CH2=CH2': -2, 'O2': -1, 'CH3CHO': 2}
# 酸化数：エチレン C -2，アセトアルデヒドのカルボニル C +1（メチル C は -3 のまま）
from fractions import Fraction
# (5)
n = 1000/44; v = n/2*22.4; print(round(v,2)); assert sig(v,3) == 255.0
print('V =', v, '→ 有効数字3桁', sig(v,3))
print('ALL OK')
