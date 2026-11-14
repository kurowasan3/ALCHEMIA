# THE-GAS-055 検証：凝縮を含む圧縮因子 Z の変化
R = 8.31e3          # Pa·L/(K·mol)
T = 363.0
RT = R*T
pv = {90: 0.70e5, 100: 1.01e5, 110: 1.43e5}
p_obs = 0.70e5
Tsel = [t for t, p in pv.items() if abs(p - p_obs) < 1e3]
print('T =', Tsel); assert Tsel == [90]
Zg = 0.993
Vl = 0.0187         # L/mol
Zl_at = lambda p: p*Vl/RT
zl = Zl_at(p_obs)
print('Z(liquid at p\') =', zl)
# (2) y = 0.50
Z_half = 0.5*Zg + 0.5*zl
print('Z(y=0.5) =', Z_half); assert round(Z_half, 3) == 0.497
# 体積からの検算：V = (1-y) Zg RT/p + y Vl
V_half = 0.5*Zg*RT/p_obs + 0.5*Vl
assert abs(p_obs*V_half/RT - Z_half) < 1e-12
# (3) V = 10.0 L
Z10 = p_obs*10.0/RT
y = (Zg - Z10)/(Zg - zl)
print('Z(10 L) =', Z10, ' y =', y); assert round(y, 2) == 0.77
# (4) 全部液体，p = 2.0e5
Z2 = Zl_at(2.0e5)
print('Z(liquid, 2.0e5) =', Z2); assert round(Z2*1e3, 1) == 1.2
# 液体の Z は p に比例（V 一定）
assert abs(Zl_at(1.0e5)*2 - Z2) < 1e-15
# 気体の体積（凝縮直前）
Vg = Zg*RT/p_obs
print('V_gas at onset =', Vg, 'L')
# (5) N2 1.00 mol を共存：水のモル分率 1/2 → 全圧 2p' で凝縮開始（理想気体近似）
p_start = p_obs/0.5
print('onset total p =', p_start); assert abs(p_start - 1.4e5) < 1
# 凝縮が進むと気相の水の分圧は p' のまま，N2 の分圧が上がり全圧は増加し続ける
for nw_gas in (1.0, 0.8, 0.5):
    Vgas = nw_gas*RT/p_obs   # 水蒸気の分圧が p' となる気相の体積
    p_tot = p_obs + 1.0*RT/Vgas
    print(' water gas', nw_gas, 'V', round(Vgas, 2), 'p_tot', round(p_tot))
print('OK')
