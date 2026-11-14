# 検証スクリプト（共通ヘッダ＋問題固有の検証）  python3 -I <this file>
import json, os, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
ID = os.path.splitext(os.path.basename(__file__))[0]
ITEM = json.load(open(f'/home/claude/chem/out2/items/{ID}/item.json', encoding='utf-8'))
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
# 小問数と答えの番号の一致
_nq = ITEM['problem'].count('\\item')
_na = len(re.findall(r'(?:^|\s)\((\d+)\)', ITEM['answer']))
assert _nq == _na or (_nq == 0 and _na == 0), (_nq, _na)
print(f'{ID}: 小問 {_nq} 個と答え {_na} 個が対応')
M = dict(H=1.0, C=12.0, N=14.0, O=16.0, Na=23.0, Mg=24.0, Al=27.0, Si=28.0, P=31.0, S=32.0, Cl=35.5, K=39.0,
         Ca=40.0, Ti=48.0, Cr=52.0, Mn=55.0, Fe=56.0, Co=59.0, Cu=63.5, Zn=65.4, Br=80.0, Sr=88.0, Ag=108.0, I=127.0, Ba=137.0, Pt=195.0)
def mass(f, masses=None):
    m = masses or M
    c, _ = parse_formula(f)
    return sum(m[e] * float(n) for e, n in c.items())
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
# ---- 問題固有の検証 ----
# 八面体錯体の立体異性体の列挙（二座配位子 en はシス位の 2 か所を占める）
from itertools import permutations, product
V = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
def mats():
    Ms = []
    for perm in permutations(range(3)):
        for sg in product((1,-1), repeat=3):
            M = [[0]*3 for _ in range(3)]
            for i in range(3): M[i][perm[i]] = sg[i]
            # 行列式
            d = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
            Ms.append((M, d))
    return Ms
MS = mats()
def act(M, i):
    v = V[i]; w = tuple(sum(M[r][c]*v[c] for c in range(3)) for r in range(3)); return V.index(w)
PERMS_ROT = [[act(M, i) for i in range(6)] for M, d in MS if d == 1]
PERMS_ALL = [[act(M, i) for i in range(6)] for M, d in MS]
def trans(i, j): return all(a == -b for a, b in zip(V[i], V[j]))
def key(conf, p):
    lab, pairs = conf
    nl = [None]*6
    for i in range(6): nl[p[i]] = lab[i]
    np_ = frozenset(frozenset((p[a], p[b])) for a, b in pairs)
    return (tuple(nl), tuple(sorted(tuple(sorted(x)) for x in np_)))
def canon(conf, group): return min(key(conf, p) for p in group)
def enumerate_isomers(mono, n_en):
    """mono: 単座配位子のリスト，n_en: en の数 → (立体異性体のリスト, 各々がキラルか)"""
    labels = mono + ['en']*(2*n_en)
    seen = {}
    for lab in set(permutations(labels)):
        pos = [i for i in range(6) if lab[i] == 'en']
        # en の位置をシス対に分ける全ての方法
        def pairings(ps):
            if not ps: yield []; return
            a = ps[0]
            for b in ps[1:]:
                if not trans(a, b):
                    rest = [x for x in ps if x not in (a, b)]
                    for r in pairings(rest): yield [(a, b)] + r
        for pr in pairings(pos):
            conf = (lab, pr)
            c = canon(conf, PERMS_ROT)
            if c not in seen:
                seen[c] = conf
    res = []
    for c, conf in seen.items():
        chiral = canon(conf, PERMS_ALL) != canon(conf, PERMS_ALL) or None
        # キラル判定：鏡像（反転）を回転で重ねられるか
        lab, pr = conf
        inv = [V.index(tuple(-x for x in V[i])) for i in range(6)]
        mirror = key(conf, inv)
        mconf = (mirror[0], [tuple(x) for x in mirror[1]])
        chiral = canon(mconf, PERMS_ROT) != c
        res.append((conf, chiral))
    return res
def summary(mono, n_en):
    r = enumerate_isomers(mono, n_en)
    n_stereo = len(r); n_chiral = sum(1 for _, ch in r if ch)
    n_geo = n_stereo - n_chiral // 2
    return n_stereo, n_geo, n_chiral
S = summary
assert S([], 3) == (2, 1, 2)                                   # (1) [Co(en)3]3+
assert S(['Cl','Cl'], 2) == (3, 2, 2)                          # (2) [CoCl2(en)2]+：トランス 1 + シス鏡像 2
assert S(['Cl','Cl','NH3','NH3'], 1) == (4, 3, 2)              # (3)
# (4) 候補
cand = {'[CoCl(NH3)5]': S(['Cl']+['NH3']*5, 0), '[CoCl(NH3)3(en)]': S(['Cl','NH3','NH3','NH3'], 1), '[CoCl(NH3)(en)2]': S(['Cl','NH3'], 2)}
print(cand)
ok = [k for k, v in cand.items() if v[0] == 2 and v[2] == 0]; assert ok == ['[CoCl(NH3)3(en)]']
# 検算：既知の値（MA2B2C2 は 6 種類，[CoCl2(NH3)4] 2 種類）
assert S(['A','A','B','B','C','C'], 0)[0] == 6 and S(['Cl','Cl']+['NH3']*4, 0)[0] == 2
# (3) の各幾何異性体の内訳：Cl どうしトランス / NH3 どうしトランス / 両方シス（キラル）
res = enumerate_isomers(['Cl','Cl','NH3','NH3'], 1)
kinds = []
for (lab, pr), ch in res:
    cl = [i for i in range(6) if lab[i] == 'Cl']; nh = [i for i in range(6) if lab[i] == 'NH3']
    kinds.append((trans(*cl), trans(*nh), ch))
print(sorted(kinds)); assert sorted(kinds) == [(False, False, True), (False, False, True), (False, True, False), (True, False, False)]
# (4) AgCl
n = 2.87/(108+35.5); print('AgCl mol', n); assert abs(n/0.0100 - 2) < 0.01
print('ALL OK')
