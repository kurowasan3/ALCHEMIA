# POL-SUG-050 検証：アミノ化 α-シクロデキストリン（6 単位）の位置異性体と部分加水分解断片
from itertools import combinations
def canon(r): return min(r[i:] + r[:i] for i in range(len(r)))       # 回転のみ同一視（向きは区別）
def frags(r, Ls):
    n = len(r); return {''.join(r[(i+k) % n] for k in range(L)) for L in Ls for i in range(n)}
def isomers(n, k): return sorted({canon(''.join('N' if i in c else 'G' for i in range(n))) for c in combinations(range(n), k)})
res = {k: isomers(6, k) for k in (1, 2, 3)}
print(res); assert [len(res[k]) for k in (1, 2, 3)] == [1, 3, 4]
# 向きを区別しない（裏返しも同一視）ときの数
canon2 = lambda r: min(canon(r), canon(r[::-1]))
print('裏返しも同一視:', len({canon2(s) for s in res[3]})); assert len({canon2(s) for s in res[3]}) == 3
assert canon('NNGNGG'[::-1]) == canon('NNGGNG')
for s in res[3]:
    print(s, '単糖〜三糖', len(frags(s, (1, 2, 3))), ' 単糖〜四糖', len(frags(s, (1, 2, 3, 4))), sorted(frags(s, (4,))))
cnt = {s: len(frags(s, (1, 2, 3))) for s in res[3]}
assert cnt == {'GGGNNN': 12, 'GGNGNN': 12, 'GGNNGN': 12, 'GNGNGN': 6}
X = [s for s in res[3] if 'NNN' not in frags(s, (3,)) and 'NNGN' in frags(s, (4,))]
print('X', X); assert X == ['GGNNGN'] and canon('NNGNGG') == 'GGNNGN'
print('X の種類（単糖〜四糖）', len(frags(X[0], (1, 2, 3, 4))))
# 二置換体（参考）
print({s: len(frags(s, (1, 2, 3))) for s in res[2]})
print('OK')
