ITEM = '/home/claude/chem4/out/items/INO-PPT-036/item.json'
# 検証スクリプト（グループ g07）  python3 -I この.py
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
print('小問数', n_items, '答え', n_ans)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from itertools import permutations
# 溶解性の表（高校教科書の標準的記述）: water, HNO3, gas_with_acid, NaOH(excess), NH3(excess), color_in_water
P = {
 'NaCl':   dict(w=True,  a=True,  g=False, b=True,  n=True,  blue=False),
 'Na2CO3': dict(w=True,  a=True,  g=True,  b=True,  n=True,  blue=False),
 'CuSO4':  dict(w=True,  a=True,  g=False, b=None,  n=None,  blue=True),
 'CaCO3':  dict(w=False, a=True,  g=True,  b=False, n=False, blue=False),
 'ZnO':    dict(w=False, a=True,  g=False, b=True,  n=None,  blue=False),
 'BaSO4':  dict(w=False, a=False, g=False, b=False, n=False, blue=False),
 'AgCl':   dict(w=False, a=False, g=False, b=None,  n=True,  blue=False),
}
L = 'ABCDEFG'
def ok(asg):
    p = {k: P[v] for k, v in zip(L, asg)}
    # 実験1
    if [k for k in L if p[k]['w']] != ['A','B','C']: return False
    if [k for k in L if p[k]['blue']] != ['C']: return False
    # 実験2
    if [k for k in L if not p[k]['a']] != ['F','G']: return False
    if [k for k in L if p[k]['g']] != ['B','D']: return False
    # 実験3（D,E のみ）
    if not (p['E']['b'] is True and p['D']['b'] is False): return False
    # 実験4（F,G のみ）
    if not (p['G']['n'] is True and p['F']['n'] is False): return False
    return True
sols = [a for a in permutations(P) if ok(a)]
print(sols)
assert sols == [('NaCl','Na2CO3','CuSO4','CaCO3','ZnO','BaSO4','AgCl')]
# (5) 実験1・2 だけで区別できない組
sig12 = {k: (P[k]['w'], P[k]['blue'], P[k]['a'], P[k]['g']) for k in P}
from collections import defaultdict
grp = defaultdict(list)
for k, v in sig12.items(): grp[v].append(k)
same = [g for g in grp.values() if len(g) > 1]
print('実験1,2で同じ結果:', same)
assert same == [['BaSO4', 'AgCl']]
# 図鑑
z = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8'))['results']
assert z['Zn^2+']['NaOH-ex']['f'] == '[Zn(OH)4]^2-'
assert 'アンモニア水に溶ける' in z['Ag+']['Cl']['note']
print('図鑑 OK'); print('すべて OK')
