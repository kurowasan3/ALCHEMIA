# THE-SOL-084 検証：蒸気圧降下（液面差 h）でみるオリゴ糖の加水分解
from fractions import Fraction as F
k = F(1, 100)                        # n = 0.0100 h（mol, mm）
nS = k*1; M = F(828, 100)/nS; n = (M - 18)/162
print('S', nS, 'mol  M', M, ' n', n); assert M == 828 and n == 5
bonds = (n - 1)*nS                   # 鎖状：グリコシド結合 n-1 個
cut = k*3 - k*1                      # HCl は A にも同量入れたので差に効かない
print('切断', cut, '全結合', bonds, '割合', cut/bonds); assert cut/bonds == F(1, 2)
HCl = F(5, 1000); h_Bonly = (k*3 + 2*HCl)/k; print('HCl を B だけ', h_Bonly); assert h_Bonly == 4
h_full = (n*nS)/k; print('完全加水分解', h_full); assert h_full == 5
# 環状（グルコース 6 単位，分子量 972）9.72 g：1 分子 1 か所の切断では粒子数不変
nC = F(972, 100)/972; print('環状', nC, 'h', nC/k, '1か所ずつ切断後も', nC/k)
assert nC/k == 1
print('OK')
