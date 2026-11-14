# THE-BND-054 検証：共鳴による N-O 結合の平均の多重度と結合距離の順序
from fractions import Fraction as Fr
# 八隅則を満たす電子式（共鳴構造）を，各 N-O 結合の多重度の組として列挙
res = {
 'NO2+': [(2, 2)],
 'NO2-': [(2, 1), (1, 2)],
 'NO3-': [(2, 1, 1), (1, 2, 1), (1, 1, 2)],
 'NH2OH': [(1,)],
}
avg = {k: Fr(sum(sum(s) for s in v), len(v)*len(v[0])) for k, v in res.items()}
for k, v in res.items():
    # 等価な結合の平均：各結合位置で共鳴構造の平均をとる → すべて同じになることを確認
    pos = [Fr(sum(s[i] for s in v), len(v)) for i in range(len(v[0]))]
    assert len(set(pos)) == 1
    print(k, '共鳴構造の数', len(v), ' 平均の多重度', pos[0])
# 電子数の確認（価電子の総数）
val = {'N': 5, 'O': 6, 'H': 1}
tot = {'NO2+': 5+12-1, 'NO2-': 5+12+1, 'NO3-': 5+18+1}
# 結合電子 + 非共有電子 = 総数 となるかを，代表的な電子式で確認
# NO2+: O=N=O 共有 8 電子，O の非共有 2 対ずつ 8 → 16
assert tot['NO2+'] == 8 + 8
# NO2-: O=N-O^- 共有 6 電子，N 非共有 1 対 2，O(二重) 2 対 4，O(単) 3 対 6 → 18
assert tot['NO2-'] == 6 + 2 + 4 + 6
# NO3-: 共有 8 電子，O(二重) 4，O(単)×2 12 → 24
assert tot['NO3-'] == 8 + 4 + 12
order = sorted(avg, key=lambda k: avg[k])
print('長い順（多重度の小さい順）:', order)
assert order == ['NH2OH', 'NO3-', 'NO2-', 'NO2+']
# 実測値（pm，文献の概数）との照合
exp = {'NH2OH': 146, 'NO3-': 125, 'NO2-': 124, 'NO2+': 115}
assert sorted(exp, key=lambda k: -exp[k]) == order
# HNO3：N-OH は常に単結合，残り 2 つは (2,1),(1,2) の共鳴 → 1.5
hno3 = [(1, 2, 1), (1, 1, 2)]
pos = [Fr(sum(s[i] for s in hno3), 2) for i in range(3)]
print('HNO3 の各 N-O の平均多重度', pos)
assert pos == [1, Fr(3,2), Fr(3,2)]
exp_hno3 = {'N-OH': 141, 'N-O(末端)': 121}   # 気相の実測（約 140.6, 119.9/121.1 pm）
assert exp_hno3['N-O(末端)'] < exp['NO3-'] < exp_hno3['N-OH']
print('OK')
