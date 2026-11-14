# THE-GAS-045 検証：ビクトル・マイヤー法による分子量測定（水を用いたガスビュレット）
from sympy import symbols, Eq, solve
R, T2, p, pw = 8.3e3, 300, 1.013e5, 3.6e3
m, v = 0.150, 32.0e-3
# 加熱部 A の物質量は不変 → 押し出された空気の物質量 = 試料の物質量
# ビュレット内の空気の分圧 = p - pw
M = m*R*T2/((p-pw)*v)
print('(2) M =', M)
cand = {'CH2Cl2': 12+2*1.0+2*35.5, 'CHCl3': 12+1.0+3*35.5, 'CCl4': 12+4*35.5}
best = min(cand, key=lambda k: abs(cand[k]-M))
print('    候補', cand, '→', best)
assert best == 'CHCl3' and round(M) == 119
# (3) 水蒸気圧の補正を忘れた場合
M2 = m*R*T2/(p*v)
print('(3) 補正なし M =', M2, '（小さくなる）')
assert M2 < M and round(M2) == 115
# 逆向き：CHCl3 0.150 g が水上で示すべき体積
n = 0.150/119.5
print('    期待される体積 mL =', n*R*T2/(p-pw)*1e3)
# (4) 誤差の向き（定性的に，数式で確認）
# ア：一部しか気化しない → v 小 → M 大
# イ：ビュレット内の気体の圧力が p - ρgh（< p）なのに p で計算 → v 大 → M 小
h = 1000.0   # 例：水位差 10 cm ≈ 1.0e3 Pa
v_i = n*R*T2/(p-pw-h)
assert m*R*T2/((p-pw)*v_i) < M
# エ：T1 は式に入らない → 変わらない
# オ：X 1 分子が気体 2 分子に分解 → v が 2 倍 → M 半分
assert m*R*T2/((p-pw)*2*v) < M
print('OK')
