"""ORG-DET-019 検証：図鑑 tests.json に基づく割り当ての一意性"""
import sys, json, itertools
sys.path.insert(0, '/home/claude/chem/work/g09')
from g09lib import report_equations
T = {t['id']: set(t['pos']) for t in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/tests.json'))['tests']}
cmp = ['methanol', 'propanol-2', 'acetone', 'ethyl-acetate', 'benzaldehyde', 'phenol', 'aniline', 'toluene']
miscible = {'methanol', 'propanol-2', 'acetone'}          # 水と任意の割合で混和（文献値）
solid = {'phenol'}                                         # 室温で固体（mp 41 ℃）
saponify = {'ethyl-acetate'}                               # 水に難溶で NaOH 水溶液と加熱するとけん化
def props(c):
    return dict(mis=c in miscible, solid=c in solid, na=c in T['sodium'], iod=c in T['iodoform'],
                ag=c in T['silver-mirror'], bleach=c in T['bleach'], sap=c in saponify)
sols = []
for perm in itertools.permutations(cmp):
    A, B, C, D, E, F, G, H = (props(x) for x in perm)
    P = dict(zip('ABCDEFGH', (A, B, C, D, E, F, G, H)))
    if [k for k in P if P[k]['mis']] != ['A', 'B', 'C']: continue
    if [k for k in P if P[k]['solid']] != ['D']: continue
    if [k for k in P if P[k]['na']] != ['A', 'B', 'D']: continue
    if [k for k in P if P[k]['iod']] != ['B', 'C']: continue
    if not E['ag'] or not F['bleach'] or not G['sap']: continue
    if any(H.values()): continue
    sols.append(perm)
print(sols); assert sols == [tuple(cmp[i] for i in (0, 1, 2, 5, 4, 6, 3, 7))]
# (2)：B,C の区別
B, C = 'propanol-2', 'acetone'
print('Na', B in T['sodium'], C in T['sodium'], '| iodoform', B in T['iodoform'], C in T['iodoform'])
assert (B in T['iodoform']) and (C in T['iodoform']) and (B in T['sodium']) != (C in T['sodium'])
report_equations('ORG-DET-019')
