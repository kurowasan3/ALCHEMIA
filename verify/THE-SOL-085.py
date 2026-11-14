# THE-SOL-085 検証：混合槽・2段向流・向流塔（無限段）での必要な吸収液量 W/G
import sympy as sp
f0, ff, h, A, f1, W, G = sp.symbols('f0 ff h A f1 W G', positive=True)
# (1) 混合槽1段：出口の気体と液が平衡 g=ff/h，収支 G(f0-ff)=W g
WG1 = sp.solve(sp.Eq(G*(f0-ff), W*ff/h), W)[0]/G
print('1段 W/G =', sp.simplify(WG1))
# (2) 向流塔の最小：操作線 g=(G/W)(f-ff) が平衡線 g=f/h を越えない。両端で確認
WGmin = h*(f0-ff)/f0
# 両端でのチェック：f=ff で 0 <= ff/h（常に成立），f=f0 で等号
# (3) 2段向流（各槽の出口どうしが平衡）
g2 = ff/h; g1 = f1/h
s = sp.solve([sp.Eq(G*(f1-ff), W*g2), sp.Eq(G*(f0-f1), W*(g1-g2))], [f1, f0], dict=True)[0]
print('2段：f1 =', sp.simplify(s[f1]), ' f0 =', sp.factor(s[f0]))
Av = sp.symbols('Av', positive=True)
assert sp.simplify(s[f0].subs(W, Av*h*G) - ff*(1+Av+Av**2)) == 0
vals = {f0: 0.020, ff: 0.0010, h: 2.0}
w1 = float(WG1.subs(vals)); wmin = float(WGmin.subs(vals))
Aroot = max(sp.solve(sp.Eq(1+Av+Av**2, 20), Av)); w2 = float(2.0*Aroot)
print('1段', w1, ' 2段', w2, ' 向流塔最小', wmin, ' 比 1段/塔', w1/wmin)
assert abs(w1-38) < 1e-9 and abs(wmin-1.9) < 1e-9 and abs(w2-7.775) < 1e-3
# 向流塔で W/G=1.9 のとき出口の液のモル分率 = f0/h（入口の気体と平衡）
gout = (0.020-0.0010)/1.9; print('塔出口の液 g =', gout, ' f0/h =', 0.020/2.0)
assert abs(gout-0.010) < 1e-12
# (5) 3段なら f0/ff = 1+A+A^2+A^3 = 20
A3 = [r for r in sp.solve(sp.Eq(1+Av+Av**2+Av**3, 20), Av) if r.is_real][0]; print('3段 W/G =', float(2*A3))
print('OK')
