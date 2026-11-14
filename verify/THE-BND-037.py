# 検証スクリプト THE-BND-037  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-BND-037/item.json'
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
import itertools, math
D = {'HH':436,'FF':155,'ClCl':243,'BrBr':194,'HF':567,'HCl':432,'ClF':252}
def delta(a,b,ab): return D[ab]-(D[a+a]+D[b+b])/2
dl = {'HF':delta('H','F','HF'),'HCl':delta('H','Cl','HCl'),'ClF':delta('Cl','F','ClF')}
assert dl == {'HF':271.5,'HCl':92.5,'ClF':53.0}
x = {k: 0.102*math.sqrt(v) for k,v in dl.items()}; print(x)
assert [round(x[k],2) for k in ('HF','HCl','ClF')] == [1.68,0.98,0.74]
# 丸めの境界でないこと
assert all(abs(x[k]*100 - round(x[k]*100)) < 0.45 for k in x)
# (2) 一直線上の並びの全探索：両端の差 ≒ 他の2つの和 になる並び
r = {frozenset(['H','F']):x['HF'], frozenset(['H','Cl']):x['HCl'], frozenset(['Cl','F']):x['ClF']}
ok = []
for p in itertools.permutations(['H','Cl','F']):
    ends = r[frozenset([p[0],p[2]])]; mids = r[frozenset([p[0],p[1]])]+r[frozenset([p[1],p[2]])]
    if abs(ends-mids) < 0.1: ok.append(p)
print(ok); assert set(ok) == {('H','Cl','F'),('F','Cl','H')}
# H–F で H が正 → F > H。よって F > Cl > H
import json
E = {e['sym']: e for e in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json', encoding='utf-8'))['elements']}
assert E['F']['en'] > E['Cl']['en'] > E['H']['en']
assert abs((E['Br']['en']-E['H']['en']) - 0.76) < 1e-9
# (3)
dHBr = (0.76/0.102)**2; DHBr = (436+194)/2 + dHBr; print(dHBr, DHBr)
assert sig(DHBr,3) == 371
# (4)
assert sorted(x, key=x.get, reverse=True) == ['HF','HCl','ClF']
assert abs(x['ClF']/x['HCl'] - 0.75) < 0.02

print('ALL OK')
