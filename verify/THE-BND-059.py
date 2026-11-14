# THE-BND-059 検証：形式電荷のモデルで原子の並び方と電子式を選ぶ
from itertools import product
VAL = {'H': 1, 'C': 4, 'N': 5, 'O': 6}
EN = {'H': 2.2, 'C': 2.6, 'N': 3.0, 'O': 3.4}

def lewis(atoms, skeleton, charge):
    """骨格（結合の組）に対し，八隅則（H は 2）を満たす電子式をすべて列挙し，形式電荷を返す"""
    n = len(atoms); tot = sum(VAL[a] for a in atoms) - charge
    out = []
    for orders in product((1, 2, 3), repeat=len(skeleton)):
        bonded = [0]*n
        for (i, j), o in zip(skeleton, orders):
            bonded[i] += o; bonded[j] += o
        lone = [(2 if a == 'H' else 8) - 2*b for a, b in zip(atoms, bonded)]
        if any(l < 0 for l in lone): continue
        if any(a == 'H' and b != 1 for a, b in zip(atoms, bonded)): continue
        if 2*sum(orders) + sum(lone) != tot: continue
        fc = [VAL[a] - l - b for a, l, b in zip(atoms, lone, bonded)]
        out.append((orders, fc, sum(abs(x) for x in fc)))
    return out

def best(atoms, skeleton, charge):
    L = lewis(atoms, skeleton, charge)
    m = min(x[2] for x in L)
    return m, [x for x in L if x[2] == m]

def show(name, atoms, sk, ch):
    m, B = best(atoms, sk, ch)
    print(name, '最小の|形式電荷|の和 =', m, [(o, f) for o, f, _ in B])
    return m, B

# (1) N2O：NNO と NON
m1, B1 = show('NNO', ['N', 'N', 'O'], [(0, 1), (1, 2)], 0)
m2, B2 = show('NON', ['N', 'O', 'N'], [(0, 1), (1, 2)], 0)
assert m1 == 2 and m2 == 4
# 同点の 2 つのうち，負の形式電荷が O にあるもの（N≡N-O）を優先
pref = [b for b in B1 if b[1][2] == -1]
assert pref and pref[0][0] == (3, 1)
# (2) HCN と HNC
m3, _ = show('HCN', ['H', 'C', 'N'], [(0, 1), (1, 2)], 0)
m4, _ = show('HNC', ['H', 'N', 'C'], [(0, 1), (1, 2)], 0)
assert m3 == 0 and m4 == 2
# (3) OCN- と CNO-
m5, B5 = show('OCN-', ['O', 'C', 'N'], [(0, 1), (1, 2)], -1)
m6, B6 = show('CNO-', ['C', 'N', 'O'], [(0, 1), (1, 2)], -1)
assert m5 == 1 and m6 == 3
assert any(b[1] == [-1, 0, 0] for b in B5)       # -O-C≡N（負電荷が最も電気陰性度の大きい O に）
# (4) CO：八隅則を満たすのは C≡O だけで，C -1，O +1
L = lewis(['C', 'O'], [(0, 1)], 0)
print('CO', L)
assert L == [((3,), [-1, 1], 2)]
assert EN['O'] > EN['C']     # 電気陰性度による分極は C(δ+) O(δ-)，形式電荷はその逆向き
# (5) O3：O-O-O
m7, B7 = show('O3', ['O', 'O', 'O'], [(0, 1), (1, 2)], 0)
assert m7 == 2 and all(b[1][1] == 1 for b in B7)
# 中心 O の電子対の領域：非共有電子対 1 + 結合 2 方向 = 3 → 折れ線形
for o, f, _ in B7:
    bonded_c = sum(o); lone_c = 8 - 2*bonded_c
    assert lone_c // 2 == 1
print('参考：実測の双極子モーメント N2O 0.16 D（NNO，直線），CO 0.11 D（C 側が負），O3 0.53 D（折れ線，結合角 117°）')
print('OK')
