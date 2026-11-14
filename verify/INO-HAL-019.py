# 検証スクリプト INO-HAL-019   python3 -I この.py
ITEM = '/home/claude/chemk/out/items/INO-HAL-019/item.json'
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
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
# 段階の式を 1 単位（n=1）で確認
for e in ['NH3BH3 -> NH2BH2 + H2', 'NH2BH2 -> NHBH + H2', 'NHBH -> BN + H2']:
    assert check_equation(e) is None, e
M = 14 + 10.8 + 6*1.0
nH2 = 4.00e3/2.0
V = nH2*22.4/1e3; print('(2)', V); assert sig(V,3) == 44.8
m2 = nH2/2*M/1e3; m3 = nH2/3*M/1e3; print('(3)', m2, m3); assert sig(m2,3) == 30.8 and sig(m3,3) == 20.5
pMg = 2.0/(24.3+2.0)*100; pAB = 4.0/M*100; print('(4)', pMg, pAB, 4.00/(pMg/100))
assert sig(pMg,3) == 7.60 and sig(pAB,3) == 13.0 and sig(4.00/(pMg/100),3) == 52.6 and pAB > pMg
# エタンと分子量が近いこと
assert abs(M - 30.0) < 1.0

print('検証OK')
