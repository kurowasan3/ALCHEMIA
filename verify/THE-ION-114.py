# THE-ION-114 検証：ハロゲン化物イオンの分別沈殿の完全さ（量の比だけで決まること）
KCl, KBr, KI = 1.8e-10, 5.4e-13, 8.1e-17
V0 = 0.0250; cBr = 4.00e-3; cCl = 4.00e-2
nBr = V0*cBr; nCl = V0*cCl
r = KBr/KCl; print('[Br-]/[Cl-] at AgCl onset =', r)            # 3.0e-3
for dil in (1, 10):                                              # 10 倍希釈しても同じ
    V = V0*dil + 0.010                                           # 滴下した液量に依らない
    Clc = nCl/V; Brc = Clc*r; frac_rem = Brc*V/nBr
    print('dilution', dil, 'precipitated %', 100*(1-frac_rem))
    assert abs(100*(1-frac_rem) - 97.0) < 1e-9
# 99.0 % 以上沈殿させる n_Cl/n_Br の上限
rmax_Br = 0.010/(KBr/KCl); rmax_I = 0.010/(KI/KCl)
print('max nCl/nBr =', rmax_Br, ' max nCl/nI =', rmax_I)          # 3.3, 2.2e4
assert f'{rmax_Br:.2g}' == '3.3' and f'{rmax_I:.2g}' == '2.2e+04'
# Cl- の半分が沈殿した時点で溶液中に残る Br-（[Br-]/[Cl-] は一定）
nBr_half = r*(nCl/2); print('Br- remaining (mol) =', nBr_half)   # 1.5e-6
assert f'{nBr_half:.2g}' == '1.5e-06'
# AgCl 沈殿開始時の Ag+ が AgBr の飽和を超えている（矛盾がない）ことの確認
for Vt in (0.025, 0.25):
    Ag = KCl/(nCl/(Vt)); print('Ag+ at onset', Ag)
print('OK')
