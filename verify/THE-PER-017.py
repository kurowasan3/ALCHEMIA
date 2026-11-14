# 検証スクリプト THE-PER-017  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-PER-017/item.json'
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
syms = ['Li','Be','B','C','N','O','F','Ne']
Q = [5.4, 9.3, 8.3, 11.3, 14.5, 13.6, 17.4, 21.6]
for s, q in zip(syms, Q):
    ev = E[s]['ie1']/96.485
    assert abs(ev - q) < 0.06, (s, ev, q)
R = [145,105,85,70,65,60,50]
for s, r in zip(syms, R):
    assert E[s]['radius'] == r, (s, E[s]['radius'])
assert all(R[i] > R[i+1] for i in range(6))
P = [1,2,3,4,5,6,7,0]
# 価電子数：族から
val = lambda g: 0 if g == 18 else (g if g <= 2 else g-10)
assert P == [val(E[s]['group']) for s in syms]
# Q の落ち込み位置
drops = [i+3 for i in range(7) if Q[i+1] < Q[i]]
assert drops == [4, 7], drops
S = [181,1287,2076,3500,-210,-218,-220,-249]
assert S.index(max(S)) == 3
# Na vs Li
assert val(E['Na']['group']) == val(E['Li']['group'])
assert E['Na']['ie1'] < E['Li']['ie1'] and E['Na']['radius'] > E['Li']['radius']
assert 98 < 181
print('Na eV', round(E['Na']['ie1']/96.485,2), 'radius', E['Na']['radius'])

print('ALL OK')
