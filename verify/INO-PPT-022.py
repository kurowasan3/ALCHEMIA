# 検証スクリプト（自動生成ヘッダ＋問題固有の検証）  python3 -I <this file>
import json, os, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
ID = os.path.splitext(os.path.basename(__file__))[0]
ITEM = json.load(open(f'/home/claude/chem/out/items/{ID}/item.json', encoding='utf-8'))
def _ce(s):
    out = []; i = 0
    while (i := s.find('\\ce{', i)) >= 0:
        j = i + 4; dep = 1
        while dep:
            dep += {'{': 1, '}': -1}.get(s[j], 0); j += 1
        out.append(s[i+4:j-1]); i = j
    return out
EQS = []
for _k in ('problem', 'answer', 'solution'):
    for _c in _ce(ITEM[_k]):
        if len(re.findall(r'<=>|->|<-', _c)) == 1 and not re.search(r'\\boxed|\\phantom|\\text|\bR\b|\bn\b', _c):
            EQS.append(_c)
for _c in EQS:
    _r = check_equation(_c)
    assert _r is None, (_c, _r)
print(f'{ID}: 反応式 {len(EQS)} 本はすべて原子数・電荷がつり合っている')
M = dict(H=1.0, He=4.0, C=12.0, N=14.0, O=16.0, Na=23.0, Mg=24.0, Al=27.0, Si=28.0, P=31.0, S=32.0, Cl=35.5, K=39.0,
         Ca=40.0, Ti=48.0, Cr=52.0, Mn=55.0, Fe=56.0, Co=59.0, Cu=63.5, Zn=65.4, Ga=69.7, Br=80.0, Ag=108.0, I=127.0, Pt=195.0, Au=197.0, Pb=207.0)
def mass(f, masses=None):
    m = masses or M
    c, _ = parse_formula(f)
    return sum(m[e] * float(n) for e, n in c.items())
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
# ---- 問題固有の検証 ----
# 頂点の置換群で配置の同値類を数える（鏡像は区別しない＝回転＋鏡映）
import itertools
import numpy as np
def sym_group(P):
    P = np.array(P, float); n = len(P); perms = set()
    base = next(b for b in itertools.combinations(range(n), 2) if np.linalg.norm(np.cross(P[b[0]], P[b[1]])) > 1e-6)
    B = np.array([P[base[0]], P[base[1]], np.cross(P[base[0]], P[base[1]])])
    for img in itertools.permutations(range(n), 2):
        for sgn in (1, -1):
            I = np.array([P[img[0]], P[img[1]], sgn * np.cross(P[img[0]], P[img[1]])])
            try: T = np.linalg.solve(B, I)
            except np.linalg.LinAlgError: continue
            if not np.allclose(T @ T.T, np.eye(3), atol=1e-6): continue
            perm = []
            for q in P @ T:
                d = np.linalg.norm(P - q, axis=1); j = int(np.argmin(d))
                if d[j] > 1e-6: break
                perm.append(j)
            else:
                if len(set(perm)) == n: perms.add(tuple(perm))
    return perms
def count(P, labels):
    G = sym_group(P); seen = set(); c = 0
    for arr in set(itertools.permutations(labels)):
        if arr in seen: continue
        c += 1
        for g in G: seen.add(tuple(arr[g[i]] for i in range(len(arr))))
    return c, len(G)
hexagon = [(np.cos(k*np.pi/3), np.sin(k*np.pi/3), 0) for k in range(6)]
prism = [(np.cos(k*2*np.pi/3), np.sin(k*2*np.pi/3), z) for z in (0.8, -0.8) for k in range(3)]
octa = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
square = [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
tetra = [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
res = {}
for name, P in [('hex', hexagon), ('prism', prism), ('octa', octa)]:
    for lab in ('AAAABB', 'AAABBB'):
        res[(name, lab)] = count(P, lab)
res[('square', 'AABB')] = count(square, 'AABB'); res[('tetra', 'AABB')] = count(tetra, 'AABB')
for k, v in res.items(): print(k, '異性体', v[0], '群の位数', v[1])
assert [res[(s, 'AAAABB')][0] for s in ('hex', 'prism', 'octa')] == [3, 3, 2]
assert [res[(s, 'AAABBB')][0] for s in ('hex', 'prism', 'octa')] == [3, 3, 2]
assert res[('square', 'AABB')][0] == 2 and res[('tetra', 'AABB')][0] == 1
assert [res[(s, 'AAAABB')][1] for s in ('hex', 'prism', 'octa')] == [12, 12, 48]
print('ALL OK')
