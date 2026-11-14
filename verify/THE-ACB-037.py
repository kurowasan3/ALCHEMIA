# THE-ACB-037 検証：濃い塩酸の pH と NaOH 滴下の上限
import math
log2 = 0.30
pH0 = -math.log10(2.0); print('pH0', pH0)
H = (2.0*10.0 - 0.50*30.0)/(40.0); print('[H+] after 30 mL', H, 'pH', -math.log10(H))
assert abs(H - 0.125) < 1e-12 and round(-math.log10(H), 2) == 0.90
Veq = 2.0*10.0/0.50; print('Veq', Veq)
lim = 14 + math.log10(0.50); print('upper limit pH', lim)
assert round(lim, 2) == 13.70
# 大過剰（V mL 加えた）でも [OH-] = (0.5V - 20)/(10+V) < 0.5
for V in [100, 1000, 1e5]:
    print(V, 14 + math.log10((0.5*V - 20)/(10 + V)))
lim1 = 14 + math.log10(1.0); print('with 1.0 M NaOH limit', lim1, 'Veq', 20.0/1.0)
print('OK')
