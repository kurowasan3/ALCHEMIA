# 検証スクリプト THE-ACB-028  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/THE-ACB-028/item.json'
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
def react(st):
    H, C, B, co2 = st
    x = min(H, C); H -= x; C -= x; B += x
    y = min(H, B); H -= y; B -= y; co2 += y
    return [H, C, B, co2]
def order1(c, h, n=20000):
    st = [0, c, 0, 0]; start = None
    for i in range(n):
        st[0] += h/n; st = react(st)
        if start is None and st[3] > 1e-12: start = (i+1)/n*h
    return st, start
def order2(c, h, n=20000):
    st = [h, 0, 0, 0]; stop = None; prev = 0
    for i in range(n):
        st[1] += c/n; st = react(st)
        if stop is None and st[0] < 1e-12: stop = (i+1)/n*c
    return st, stop
s1, start = order1(0.0200, 0.0300)
s2, stop = order2(0.0200, 0.0300)
assert abs(s1[3]-0.0100) < 1e-6 and abs(s1[2]-0.0100) < 1e-6 and s1[1] < 1e-9
assert abs(s2[3]-0.0150) < 1e-6 and abs(s2[1]-0.0050) < 1e-6 and s2[2] < 1e-6
assert abs(start/0.300*1000 - 66.7) < 0.1 and abs(stop/0.200*1000 - 75.0) < 0.1, (start, stop)
# 電荷
assert abs(0.0400 - (0.0300 + s1[2])) < 1e-6 and abs(0.0400 - (0.0300 + 2*s2[1])) < 1e-6
# (4) 全探索（解析式）
def I(a,b): return 0 if b <= a else (b-a if b <= 2*a else a)
def II(a,b): return b/2 if b <= 2*a else a
sols = [(a/10000, b/10000) for a in range(1, 600) for b in range(1, 1200)
        if abs(I(a/10000,b/10000)-0.0100) < 1e-9 and abs(II(a/10000,b/10000)-0.0125) < 1e-9]
assert sols == [(0.015, 0.025)], sols
t1, _ = order1(0.0150, 0.0250); t2, _ = order2(0.0150, 0.0250)
assert abs(t1[3]-0.0100) < 1e-6 and abs(t2[3]-0.0125) < 1e-5
for e in ['Na2CO3 + HCl -> NaHCO3 + NaCl', 'NaHCO3 + HCl -> NaCl + H2O + CO2', 'Na2CO3 + 2HCl -> 2NaCl + H2O + CO2']:
    assert check_equation(e) is None
print('start mL', start/0.3*1000, 'stop mL', stop/0.2*1000)

print('ALL OK')
