# INO-OXS-020 verification: sulfur chlorides and S4N4 synthesis
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for eq in ('S2Cl2 + Cl2 -> 2SCl2', '6SCl2 + 16NH3 -> S4N4 + 2S + 12NH4Cl'):
    r = check_equation(eq); print(eq, r); assert r is None
S, Cl, N, H = 32, 35.5, 14, 1.0
wS = 2 * S / (2 * S + 2 * Cl); print('S% in S2Cl2', wS); assert round(wS * 100, 1) == 47.4
assert 120 <= 2 * S + 2 * Cl < 150
nA = 13.5 / 135; mB = 13.5 + nA * 71; print('mass B', mB); assert abs(mB - 20.6) < 0.01
# electron balance: 4 S +2->+3 lose 4 ; 2 S +2->0 gain 4
assert 4 * 1 == 2 * 2
M = 4 * S + 4 * N; n = 9.2 / M; print(M, n)
mSCl2 = n * 6 * (S + 2 * Cl); mNH3 = n * 16 * (N + 3 * H); print(mSCl2, mNH3)
assert round(mSCl2) == 31 and round(mNH3) == 14
print('OK')
