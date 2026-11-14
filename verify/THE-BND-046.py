# 検証スクリプト THE-BND-046  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-BND-046/item.json'
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
E = {e['sym']: e for e in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json', encoding='utf-8'))['elements']}
tbl = {'A':('Mg',[2,8,2]), 'B':('O',[2,6]), 'C':('Cl',[2,8,7]), 'D':('C',[2,4]), 'E':('K',[2,8,8,1])}
for k,(s,sh) in tbl.items():
    assert E[s]['shells'] == sh, (s, E[s]['shells'])
    assert sum(sh) == E[s]['z']
charge = {'Mg':2,'O':-2,'Cl':-1,'K':1}
ne = {s: E[s]['z']-charge[s] for s in charge}
assert ne == {'Mg':10,'O':10,'Cl':18,'K':18}
ionic = {'MgO':('Mg',1,'O',1), 'MgCl2':('Mg',1,'Cl',2), 'K2O':('K',2,'O',1), 'KCl':('K',1,'Cl',1)}
for f,(a,na,b,nb) in ionic.items():
    assert charge[a]*na + charge[b]*nb == 0, f
same = [f for f,(a,na,b,nb) in ionic.items() if ne[a]==ne[b]]
assert same == ['MgO','KCl'], same
# 共有結合：原子価の和
val = {'C':4,'O':2,'Cl':1}
assert val['C'] == 2*val['O'] and val['C'] == 4*val['Cl']

print('ALL OK')
