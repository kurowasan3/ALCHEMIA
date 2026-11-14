# 検証スクリプト THE-THC-070  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-THC-070/item.json'
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|\$', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式の検査', n, '本 OK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

Eg = {'X': 2.48, 'Y': 1.55, 'Z': 0.992}
edge = {k: 1240 / v for k, v in Eg.items()}
print(edge)
assert sig(edge['X'], 3) == 500 and sig(edge['Y'], 3) == 800 and sig(edge['Z'], 3) == 1250
def stack(order, lo=300, hi=2500):
    rng = {}; cur = lo   # 届いている光の最短波長
    for k in order:
        e = edge[k]
        if e > cur: rng[k] = (cur, e); cur = e
        else: rng[k] = None
    return rng, (cur, hi)
r2 = stack('XYZ'); r3 = stack('ZYX'); print(r2, r3)
assert r2[0] == {'X': (300, 500.0), 'Y': (500.0, 800.0), 'Z': (800.0, 1250.0)} and abs(r2[1][0] - 1250) < 1e-9
assert r3[0]['Y'] is None and r3[0]['X'] is None
E400 = 1240 / 400
assert abs(E400 - 3.10) < 1e-12 and sig(2.48/E400*100, 2) == 80 and sig(0.992/E400*100, 2) == 32
assert sig(1240/1550, 3) == 0.800
print('hc/eV', 6.63e-34*3.00e8/1.60e-19*1e9)

print('ALL OK')
