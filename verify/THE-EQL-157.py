# THE-EQL-157 検証：競争実験による速度論的同位体効果（分子間・分子内）
import math
from sympy import symbols, Function, dsolve, Eq, exp, solve, nsimplify
kH_over_kD = 6.0
fH = 0.250; fD = fH**(1/kH_over_kD)
consumed = 1 - (fH + fD)/2
print('fD', fD, '消費', consumed, '残りの D:H', fD/fH)
# 問題の数値から逆算
c = 0.478; ratio = 3.17
fH2 = 2*(1-c)/(1+ratio); fD2 = ratio*fH2
print('fH', fH2, 'fD', fD2, 'kH/kD', math.log10(fH2)/math.log10(fD2))
print('与えた近似で', (-0.602)/(-0.100))
# 数値積分で関係式 log fH = (kH/kD) log fD を確認（酸化剤の濃度変化は任意）
H, D, Ox = 1.0, 1.0, 1.1; kH, kD = 6.0, 1.0; dt = 1e-5
while (H + D)/2 > 1 - consumed:
    ox = Ox
    dH = -kH*H*ox*dt; dD = -kD*D*ox*dt
    H += dH; D += dD; Ox += dH + dD   # 酸化剤も消費（量的関係は任意の比で可）
print('数値積分 fH', H, 'fD', D, 'ln比', math.log(H)/math.log(D))
# 分子内競争：C6H5D のニトロ化
for name, kpos_D in (('付加が律速', 1.0), ('C-H 開裂が律速', 1/kH_over_kD)):
    keep = 5*1.0; lose = 1*kpos_D
    print(name, 'D を保持 : 失う =', keep/lose, ': 1 ，保持の割合', keep/(keep+lose))
