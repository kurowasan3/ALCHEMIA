# THE-ION-118 加水分解による多核イオンの生成と pH：独立計算
from sympy import symbols, nsolve, solve, log, N, Rational, sqrt
K1 = 1.0e-8; K22 = 4.0e-11
m, h = symbols('m h', positive=True)
# (1) pH 5.5 に保つ：全濃度 C の金属の何%が二核イオンか
H = 10**-5.5
for C in (0.10, 1.0e-3):
    sols = solve((1 + K1/H)*m + 2*K22*m**2/H**2 - C, m)
    mm = [s for s in sols if s.is_real and s > 0][0]
    frac_d = 2*K22*mm**2/H**2/C; frac_h = K1/H*mm/C
    print(f'C={C}: [M2+]={float(mm):.4g}, dimer fraction={float(frac_d):.4f}, MOH+ fraction={float(frac_h):.5f}')
# (2)(3) 緩衝なし：電荷収支 h = [MOH+] + 2[M2(OH)2] （水の電離は無視），[M2+]≈C
for C in (0.10, 1.0e-3):
    hh = nsolve(h**3 - K1*C*h - 2*K22*C**2, h, 1e-4)
    # 近似 [M2+]=C の確認
    dimerM = 2*K22*C**2/hh**2; mohp = K1*C/hh
    print(f'C={C}: h={float(hh):.4g}, pH={-float(N(log(hh,10))):.3f}; M in dimer={float(dimerM):.3g} ({float(dimerM/C):.2e}), MOH+={float(mohp):.3g}')
    # 単核のみ・二核のみの近似
    h_mono = sqrt(K1*C); h_di = (2*K22*C**2)**(Rational(1,3))
    print('   mono-only pH', -float(N(log(h_mono,10))), ' dimer-only pH', -float(N(log(h_di,10))))
# 濃度 1/100 での pH 変化
p1 = -float(N(log(nsolve(h**3 - K1*0.1*h - 2*K22*0.01, h, 1e-4),10)))
p2 = -float(N(log(nsolve(h**3 - K1*1e-3*h - 2*K22*1e-6, h, 5e-6),10)))
print('ΔpH over 2 decades =', p2-p1, ' mono model 1.00, dimer model', 4/3)
# 水の電離の無視の妥当性
print('Kw/h at C=1e-3 :', 1e-14/5.1e-6)
# (4) 二核由来の H+ が単核由来の 1/10 になる濃度（h≈sqrt(K1 C) として）
C = symbols('C', positive=True)
ratio = 2*K22*C**2/(K1*C*sqrt(K1*C))   # (2[dimer])/[MOH+]  with h = sqrt(K1 C)
Cb = solve(ratio - Rational(1,10), C)[0]; print('C boundary =', float(Cb), ' h there =', float(sqrt(K1*Cb)))
# 反復解（問題文のヒントどおり）
for C0 in (0.1, 1e-3):
    x = (2*K22*C0**2)**(1/3)
    for i in range(6): x = (K1*C0*x + 2*K22*C0**2)**(1/3)
    print('iter', C0, x)
