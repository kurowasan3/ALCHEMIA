# INO-CSI-013 検証：ケイ素結晶の原子密度と微量元素の添加（p 型・n 型）
import sympy as sp
R = sp.Rational
a = R(543, 1000)*sp.Integer(10)**-7   # cm
n_cell = 8*R(1, 8) + 6*R(1, 2) + 4
assert n_cell == 8
nSi = n_cell/a**3; print('Si atoms per cm3', float(nSi)); assert f"{float(nSi):.1e}" == '5.0e+22'
holes = R(25, 10)*sp.Integer(10)**17
ratio = holes/R(50, 10)/sp.Integer(10)**22   # 有効数字に合わせ 5.0e22 を使用
print('B:Si', float(ratio)); assert float(ratio) == 5.0e-6
mB = R(100, 100)/R(281, 10)*ratio*R(108, 10)
print('mass B in 1.00 g Si (g)', float(mB)); assert f"{float(mB):.1e}" == '1.9e-06'
groups = {'Ga': 13, 'P': 15, 'Sb': 15, 'In': 13, 'Al': 13}
print('n-type dopants', [k for k, g in groups.items() if g == 15])
