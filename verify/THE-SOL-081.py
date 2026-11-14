# THE-SOL-081 検証：2成分混合物の凝固点降下と独立な情報
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import symbols, Eq, solve, Matrix
kf = 1.85; W = 0.200
dT = round(kf*(1.80/180 + 3.20/92)/W, 3); print('dT (given)', dT)
assert dT == 0.414
n_tot = dT/kf*W
n_glc = 2*(0.715/143)          # 2.50 g 分の Cu2O → 5.00 g 分
MB = (5.00 - n_glc*180)/(n_tot - n_glc)
print('glucose frac', n_glc*180/5.00, 'MB', MB)
assert abs(n_glc*180/5.00 - 0.36) < 1e-9 and round(MB) == 92
# 独立性：未知数 f（グルコースの質量分率），u=1/MB。ΔT 測定は W_X*(f/180 + (1-f)u) の形
f, u, w1, w2, wg = symbols('f u w1 w2 wg', positive=True)
g1 = w1*(f/180 + (1-f)*u)               # 質量 w1 の X
g2 = w2*(f/180 + (1-f)*u)               # 質量 w2 の X（①）
g3 = w1*(f/180 + (1-f)*u) + wg/180      # X＋既知のグルコース（②）
J = lambda a, b: Matrix([a, b]).jacobian([f, u]).det().simplify()
print('det ①', J(g1, g2), ' det ②', J(g1, g3))
assert J(g1, g2) == 0 and J(g1, g3) == 0
print('OK')
