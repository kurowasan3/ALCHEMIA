# THE-EQL-131 検証：生成物の一部が凝縮する定積一次分解と全圧
import math
Pv = 0.30  # C の飽和蒸気圧（P0 単位）
def P_total(x, P0=1.0):
    pC = x*P0
    if pC <= Pv: return P0*(1 + x)
    return P0*(1 - x) + x*P0 + Pv
# 10 分で 1.20 P0 → x = 0.20（凝縮前：C の分圧 0.20 < 0.30）
x10 = 0.20; assert abs(P_total(x10) - 1.20) < 1e-12
r = 0.80  # 10 分ごとに残る割合
log2, log7 = 0.301, 0.845
log08 = 3*log2 - 1; log07 = log7 - 1
t_c = 10*log07/log08
print('condensation starts at', t_c, 'exact', 10*math.log(0.7)/math.log(0.8))
assert round(t_c) == 16
for t in [0, 10, 16, 20, 30, 60]:
    x = 1 - r**(t/10); print(t, round(x, 4), round(P_total(x), 4))
x30 = 1 - r**3
print('30 min: remaining A', r**3, 'liquid C (per A0)', x30 - Pv)
assert abs(r**3 - 0.512) < 1e-12 and abs(x30 - Pv - 0.188) < 1e-12
# 30 分の分圧：A 0.512, B 0.488, C 0.30
print('pA,pB,pC', r**3, x30, Pv, 'sum', r**3 + x30 + Pv)
# 初圧の条件：凝縮しない ⇔ P0' ≤ Pv（全部分解しても C の分圧 ≤ Pv）
P0p = 0.20
print("P0'=0.20, 30 min total", P_total(x30, P0p) if x30*P0p <= Pv else None, P0p*(1 + x30))
assert abs(P0p*(1 + x30) - 0.2976) < 1e-12
# 初圧 0.50 P0 での凝縮の開始：x*0.50 = 0.30 → x = 0.60 → 残り 0.40 → t = 10 log0.4/log0.8
t2 = 10*math.log(0.4)/math.log(0.8); t2g = 10*(2*log2 - 1)/log08
print('P0=0.5: start', t2, t2g)
assert round(t2g) == 41
print('final total pressure (t->inf) for P0=1:', 1 + Pv, 'for 0.5 P0:', 0.5 + Pv)
print('OK')
