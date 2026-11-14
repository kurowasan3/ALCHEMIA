# INO-AMP-019 検証：アルミニウムの陽極酸化（電流効率・酸化皮膜の溶解率・ピリング-ベッドワース比）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Eq, Rational
Al, O, F = 27.0, 16.0, 9.65e4
Al2O3 = 2*Al + 3*O
I, t = 20.0e-3, 1930
Q = I*t; ne = Q/F
print('Q =', Q, 'C  n(e-) =', ne)
# 正方向の設定（真の値）：f=0.900, d=0.200 → 測定値をつくる
f0, d0 = 0.900, 0.200
nAl0 = f0*ne/3
loss0 = nAl0*Al*1e3                              # mg（皮膜をはがした後の質量減少）
gain0 = (-nAl0*Al + (1-d0)*nAl0/2*Al2O3)*1e3     # mg（電解後の質量増加）
print('測定値：はがした後の減少 %.3f mg，電解後の増加 %.3f mg' % (loss0, gain0))
assert round(loss0,2) == 3.24 and round(gain0,2) == 1.66
# 逆に問題文の値から求める
loss, gain = 3.24e-3, 1.66e-3
nAl = loss/Al
f = 3*nAl/ne
mox = gain + loss
nox = mox/Al2O3
d = 1 - nox/(nAl/2)
print('(2) n(Al) =', nAl, ' 効率 f =', f, ' O2 =', (1-f)*ne/4, 'mol')
print('(3) 残った Al2O3 =', mox*1e3, 'mg =', nox, 'mol  溶けた割合 d =', d)
assert abs(f-0.900) < 0.001 and abs(d-0.20) < 0.005
# 記号でも確認
fs, ds = symbols('f d')
sol = solve([Eq(fs*ne/3*Al, loss), Eq(-fs*ne/3*Al + (1-ds)*fs*ne/6*Al2O3, gain)], [fs, ds])
print('   連立の解', sol)
# (4) PBR
pbr_al = (Al2O3/3.95)/(2*Al/2.70)
pbr_mg = (40.3/3.58)/(24.3/1.74)
print('(4) PBR Al =', pbr_al, ' Mg =', pbr_mg)
assert round(pbr_al,2) == 1.29 and round(pbr_mg,2) == 0.81
for e in ['2Al + 3H2O -> Al2O3 + 6H+ + 6e-', '2H2O -> O2 + 4H+ + 4e-', 'Al2O3 + 6H+ -> 2Al^3+ + 3H2O']:
    assert check_equation(e) is None, e; print('式OK', e)
nn = 1.66e-3/24.0
print('誤った考え：n =', nn, ' 効率', 3*nn/ne)
assert round(3*nn/ne,2) == 0.52
print('OK')
