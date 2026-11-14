# THE-SOL-093 検証：ベンゼン中で完全に二量体 → 分子量，水中の凝固点降下 → 電離度・電離定数，希釈の効果
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
# 電離の式は一般式（HA）なので chemlib の対象外
from sympy import Rational as R, symbols, solve, N
Kb, Kw_f = R('5.12'), R('1.85')
dT1, w1, W1 = R('1.024'), R('1.89'), R('50.0')
m_dimer = dT1/Kb; n_dimer = m_dimer*W1/1000; M = w1/(2*n_dimer)
print('M =', N(M)); assert M == R('94.5')
w2, W2, dT2 = R('0.945'), R('100'), R('0.207')
m = w2/M*1000/W2; i = dT2/(Kw_f*m); alpha = i - 1
print('m =', N(m), 'alpha =', N(alpha)); assert m == R('0.1'); assert round(float(alpha), 2) == 0.12
Ka = m*alpha**2/(1 - alpha); print('Ka =', N(Ka)); assert round(float(Ka)*1e3, 1) == 1.6
# 1/10 に薄めたとき（Ka は (3) の値 1.6e-3 を使う）
a = symbols('a', positive=True)
a2 = [x for x in solve(R('0.01')*a**2/(1 - a) - R('1.6e-3'), a) if 0 < x < 1][0]
dT3 = Kw_f*R('0.01')*(1 + a2)
print('alpha(0.01) =', N(a2), 'dT =', N(dT3), 'vs dT2/10 =', N(dT2/10))
assert dT3 > dT2/10 and round(float(a2), 2) == 0.33
print('OK')
