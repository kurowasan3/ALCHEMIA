# 検証スクリプト POL-SUG-025  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/POL-SUG-025/item.json'
import sys, json, re, random
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import formula_mass
it = json.load(open(ITEM, encoding='utf-8'))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 5, (n_items, n_ans)
eqs = [b for k in ('problem', 'answer', 'solution') for b in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[k])
       if len(re.findall(r'->|<=>|<-', b)) == 1]
assert eqs == []
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
M = {'H': 1.0, 'C': 12, 'O': 16}
mt, mr, md = (formula_mass(f, M) for f in ('C10H20O6', 'C9H18O6', 'C8H16O6')); assert (mt, mr, md) == (236, 222, 208)
T, Tr, D = 0.472/mt, 3.69/mr, 0.291/md; print(T, Tr, D, T + Tr + D)
assert abs((T + Tr + D) - 3.24/162) / (3.24/162) < 0.002
nt = T - D; print('樹木状', nt); assert sig(nt, 2) == 6.0e-4
n = 2.49e4*0.100/(8.3e3*300); print('全分子', n); assert sig(n, 3) == 1.00e-3
nc = n - nt; print('環状', nc); assert sig(nc, 2) == 4.0e-4 and round(nc/n*100) == 40
assert sig(3.24/6.0e-4, 2) == 5400 and sig(3.24/1.00e-3, 2) == 3200
cu = 6.0e-4*144; print('Cu2O', cu); assert sig(cu, 2) == 0.086
# グルカンのモデルを作って e=b+1（樹木状），e=b（環状）を確かめる
# 単位 i の 1 位が単位 parent[i] の O-4 または O-6 に結合する．O-4 には高々 1 本，O-6 にも高々 1 本．
def tree(N, rng):
    occ = {0: set()}; links = []
    for i in range(1, N):
        while True:
            p = rng.randrange(i); site = rng.choice((4, 6))
            if site not in occ[p] and not (site == 6 and 4 not in occ[p] and rng.random() < 0.5):
                break
        occ[p].add(site); occ[i] = set(); links.append((i, p, site))
    return occ, len(links), N - 1   # 単位 0 が還元末端（1 位が空いている）
def ring(N, ringsize, rng):
    occ = {i: set() for i in range(N)}; nlinks = 0
    for i in range(ringsize):
        occ[(i + 1) % ringsize].add(4); nlinks += 1
    for i in range(ringsize, N):
        while True:
            p = rng.randrange(i); site = rng.choice((4, 6))
            if site not in occ[p]: break
        occ[p].add(site); nlinks += 1
    return occ, nlinks, N
rng = random.Random(1)
tot = dict(T=0, Tr=0, D=0); ntree = 0
for trial in range(400):
    if trial % 2 == 0:
        N = rng.randrange(2, 60); occ, nl, expect = tree(N, rng); kind = 'tree'; ntree += 1
    else:
        rs = rng.randrange(6, 12); N = rs + rng.randrange(0, 40); occ, nl, expect = ring(N, rs, rng); kind = 'ring'
    assert nl == expect
    e = sum(1 for s in occ.values() if len(s) == 0); b = sum(1 for s in occ.values() if len(s) == 2)
    assert sum(len(s) for s in occ.values()) == nl
    assert (e == b + 1) if kind == 'tree' else (e == b), (kind, e, b)
    tot['T'] += e; tot['D'] += b; tot['Tr'] += sum(1 for s in occ.values() if len(s) == 1)
assert tot['T'] - tot['D'] == ntree
print('モデル検証: テトラ－ジ =', tot['T'] - tot['D'], '= 樹木状の分子数', ntree)
print('全検証OK:', ITEM)
