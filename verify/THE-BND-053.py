# THE-BND-053 検証：金属・半導体の電気伝導性の温度変化と 14 族の順序
# 標準的な物性値（教科書・便覧の概数）
gap_eV = {'ダイヤモンド': 5.5, 'ケイ素': 1.12, 'ゲルマニウム': 0.67}
sigma_RT = {'ダイヤモンド': 1e-13, 'ケイ素': 4e-4, 'ゲルマニウム': 2.0}   # S/m（固有半導体，室温）
bond_kJ = {'ダイヤモンド': 346, 'ケイ素': 222, 'ゲルマニウム': 188}
metal_dRho_dT_positive = {'銅': True, 'ナトリウム': True}
semi = {'ケイ素', 'ゲルマニウム'}
# (1)
up = sorted(semi); down = sorted(metal_dRho_dT_positive)
assert up == ['ケイ素', 'ゲルマニウム'] or set(up) == semi
print('(1) 大きくなる:', up, ' 小さくなる:', down)
# (4) 伝導率の順＝禁制帯幅の逆順＝結合エネルギーの逆順
o1 = sorted(sigma_RT, key=lambda k: -sigma_RT[k])
o2 = sorted(gap_eV, key=lambda k: gap_eV[k])
o3 = sorted(bond_kJ, key=lambda k: bond_kJ[k])
assert o1 == o2 == o3 == ['ゲルマニウム', 'ケイ素', 'ダイヤモンド'], (o1, o2, o3)
print('(4)', ' > '.join(o1))
print('OK')
