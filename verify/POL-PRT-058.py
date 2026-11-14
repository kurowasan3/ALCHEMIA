# POL-PRT-058 検証：前駆体の切断・ジスルフィド結合・質量・結合の組合せ
from itertools import combinations
cys = [7, 19, 20, 33, 72, 85]
cut_after = [39, 60]; removed = range(40, 61)
chains = [(1, 39), (61, 90)]
chain_of = lambda r: 0 if r <= 39 else (None if r <= 60 else 1)
# (2) 選択肢：除去される断片にシステインがありそれが残る鎖と結合していれば除去できない
# (3) 質量
Mp, Mseg = 9950, 2330
mature = Mp + 2*18.0 - Mseg - 3*2.0
print('mature', mature); assert mature == 7650
# (4) 全ペアリング（完全マッチング）の列挙
def matchings(s):
    if not s: yield []; return
    a = s[0]
    for b in s[1:]:
        rest = [x for x in s if x not in (a, b)]
        for m in matchings(rest): yield [(a, b)] + m
allm = list(matchings(cys)); print('all matchings', len(allm)); assert len(allm) == 15
ok = [m for m in allm if sum(chain_of(a) != chain_of(b) for a, b in m) == 2]
print('2 inter + 1 intra', len(ok)); assert len(ok) == 12
# 断片：A5–A9, A17–A22, A31–A35, A70–A74, A83–A87（ジスルフィドを保ったまま）
frags = {'F1': range(5, 10), 'F2': range(17, 23), 'F3': range(31, 36), 'F4': range(70, 75), 'F5': range(83, 88)}
fr = lambda r: [k for k, v in frags.items() if r in v][0]
def complexes(m):
    # 断片を頂点，S–S を辺とする連結成分
    parent = {k: k for k in frags}
    def f(x):
        while parent[x] != x: x = parent[x]
        return x
    for a, b in m: parent[f(fr(a))] = f(fr(b))
    groups = {}
    for k in frags: groups.setdefault(f(k), set()).add(k)
    return sorted(sorted(g) for g in groups.values())
obs = [['F1', 'F4'], ['F2', 'F3', 'F5']]
consistent = [m for m in ok if complexes(m) == obs]
print('consistent', consistent); assert len(consistent) == 2
assert all((7, 72) in m for m in consistent)
print('OK')
