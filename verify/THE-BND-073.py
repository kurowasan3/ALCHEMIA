# THE-BND-073 検証：グラフェンの比表面積と積層数
NA=6.0e23; C=12; A_hex=5.2e-16  # cm^2（六角形1個に C 原子 2 個分）
m_hex=2*C/NA
ssa1=2*A_hex/m_hex      # cm^2/g（両面）
print('single layer SSA cm2/g', ssa1, 'm2/g', ssa1/1e4)
assert round(ssa1/1e4,-2)==2600
for n in [1,2,3,5]: print(n, ssa1/n/1e4)
n_ac=ssa1/1e4/1300; print('activated carbon avg layers', n_ac); assert round(n_ac)==2
# 黒鉛の 1.0 cm 立方体（密度 2.2 g/cm3）
ssa_cube=6/(2.2*1.0); print('cube cm2/g', ssa_cube, 'ratio', ssa1/ssa_cube)
# 面密度（片面）と A4（0.062 m2=620 cm2）を覆う単層グラフェンの質量
sd=m_hex/A_hex; print('g/cm2', sd, 'A4 mass g', sd*620)
assert round(sd*620,6)==round(4.769e-5,6) or abs(sd*620-4.8e-5)<1e-6
