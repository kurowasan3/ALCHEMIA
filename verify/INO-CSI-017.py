# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/INO-CSI-017/item.json'
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

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import itertools, math, networkx as nx
phi=(1+5**0.5)/2
base=[(0,1,3*phi),(1,2+phi,2*phi),(phi,2,2*phi+1)]
pts=set()
for b in base:
    for perm in [(0,1,2),(1,2,0),(2,0,1)]:
        v=[b[perm[0]],b[perm[1]],b[perm[2]]]
        for s in itertools.product([1,-1],repeat=3):
            pts.add(tuple(round(v[i]*s[i],6) for i in range(3)))
pts=list(pts); assert len(pts)==60
G=nx.Graph()
for i,j in itertools.combinations(range(60),2):
    if abs(math.dist(pts[i],pts[j])-2)<1e-4: G.add_edge(i,j)
assert G.number_of_nodes()==60 and G.number_of_edges()==90
assert all(d==3 for _,d in G.degree())
def in5(u,v):
    # u-v を含む長さ5の閉路があるか（v から u へ長さ4の単純道）
    for p in nx.all_simple_paths(G,v,u,cutoff=4):
        if len(p)==5: return True
    return False
e56=[e for e in G.edges if in5(*e)]; e66=[e for e in G.edges if not in5(*e)]
assert len(e56)==60 and len(e66)==30
M=nx.Graph(e66); assert all(M.degree(n)==1 for n in G.nodes)
print('C60: 辺90, [5,6]60, [6,6]30, [6,6] は完全マッチング')
E70=70*3//2; h=(2*E70-60)//6; assert (E70,h)==(105,25) and 70-E70+(12+h)==2
assert E70-60==45 and 70//2==35 and 45-35==10
print('C70 OK')

print('ALL OK')
