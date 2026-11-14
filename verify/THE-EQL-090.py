# THE-EQL-090 検証：NO2・NO・O3 の光定常状態
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['NO2 -> NO + O', 'O + O2 -> O3', 'NO + O3 -> NO2 + O2']:
    assert check_equation(e) is None
k1, k3 = 8.0e-3, 1.2e7          # s^-1, L/(mol s)
NO2, NO = 2.0e-9, 1.0e-9        # mol/L
O3 = k1*NO2/(k3*NO); print('(3) [O3] =', O3); assert f'{O3:.1e}' == '1.3e-09'
v = k1*NO2; print('   v1 = v3 =', v, ' check', k3*NO*O3)
# 3 式の和は 0（正味の変化なし）
# (4) NO を加えて [NO] が 2 倍、[NO2] は同じとみなす -> [O3] は 1/2
O3b = k1*NO2/(k3*2*NO); print('(4) [O3] with NO doubled =', O3b); assert abs(O3b/O3 - 0.5) < 1e-12
# 夜間 k1 = 0：O3 は NO と反応して減る一方
# (5) 定常状態に達するまでの目安：O3 の寿命 1/(k3[NO])
print('O3 lifetime =', 1/(k3*NO), 's')
print('OK')
