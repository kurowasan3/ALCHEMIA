# THE-ACB-039 検証：水和金属イオンの酸としての強さ（電荷・半径）と pH
import sys, math; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert check_equation('[Al(H2O)6]^3+ + H2O <=> [Al(OH)(H2O)5]^2+ + H3O+') is None
data={'Na+':(1,0.102,14.8),'Ca2+':(2,0.100,12.8),'La3+':(3,0.103,8.5),'Mg2+':(2,0.072,11.4),'Sr2+':(2,0.118,13.2),'Ba2+':(2,0.135,13.4),'Al3+':(3,0.054,5.0)}
# 半径がほぼ同じ組：電荷が大きいほど pKa が小さい
assert data['La3+'][2]<data['Ca2+'][2]<data['Na+'][2]
# 同じ電荷：半径が小さいほど pKa が小さい
g2=sorted(['Mg2+','Ca2+','Sr2+','Ba2+'],key=lambda k:data[k][1])
assert [data[k][2] for k in g2]==sorted(data[k][2] for k in g2)
# z^2/r との相関（順位）
for k,(z,r,p) in sorted(data.items(),key=lambda t:t[1][0]**2/t[1][1]): print(k, round(z*z/r,1), p)
def pH(c,pKa):
    Ka=10**-pKa
    x=(-Ka+math.sqrt(Ka*Ka+4*Ka*c))/2   # 水の電離は無視
    return -math.log10(x), x/c
for ion,pKa in (('Al3+',5.0),('La3+',8.5)):
    p,a=pH(0.10,pKa); print(ion,'pH',round(p,2),'電離度',a)
assert round(pH(0.10,5.0)[0],2)==3.00 and round(pH(0.10,8.5)[0],2)==4.75
# 水の電離の寄与（La: [H+]≈1.8e-5 ≫ 1e-7 で無視できる）
assert 10**-4.75>100*1e-7
print('OK')
