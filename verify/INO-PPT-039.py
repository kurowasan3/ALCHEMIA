# INO-PPT-039 検証：コバルト(III) アンミン錯塩の内圏・外圏の Cl と異性体
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['[Co(NH3)6]Cl3 -> [Co(NH3)6]^3+ + 3Cl-', '[CoCl(NH3)5]Cl2 -> [CoCl(NH3)5]^2+ + 2Cl-',
          '[CoCl2(NH3)4]Cl -> [CoCl2(NH3)4]+ + Cl-', 'Ag+ + Cl- -> AgCl']:
    assert check_equation(e) is None, e
AgCl = 107.9 + 35.5
n = 0.0100
for name, free, ions in [('X:CoCl3.6NH3', 3, 4), ('Y:CoCl3.5NH3', 2, 3), ('Z:CoCl3.4NH3', 1, 2), ('W:CoCl3.3NH3', 0, 1)]:
    print(name, 'AgCl =', round(n*free*AgCl, 2), 'g  ions per formula =', ions)
assert [round(n*k*AgCl, 2) for k in (3, 2, 1)] == [4.30, 2.87, 1.43]
# 八面体 MA4B2, MA3B3 の異性体の数を座標の置換で数える
import itertools
V = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
def rots():
    mats = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([1,-1], repeat=3):
            M = [[0]*3 for _ in range(3)]
            for i in range(3): M[i][perm[i]] = signs[i]
            # 行列式 +1 のみ（回転）
            import numpy as np
            if round(np.linalg.det(np.array(M))) == 1: mats.append(M)
    return mats
R = rots(); assert len(R) == 24
def apply(M, v): return tuple(sum(M[i][j]*v[j] for j in range(3)) for i in range(3))
def count(nB):
    seen = set(); classes = 0
    for pos in itertools.combinations(range(6), nB):
        key = frozenset(V[i] for i in pos)
        if key in seen: continue
        classes += 1
        for M in R: seen.add(frozenset(apply(M, v) for v in key))
    return classes
print('MA4B2 isomers =', count(2), ' MA3B3 isomers =', count(3))
assert count(2) == 2 and count(3) == 2
print('OK')
