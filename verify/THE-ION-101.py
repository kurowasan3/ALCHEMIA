# THE-ION-101 検証：炭酸水素ナトリウム水溶液を大気と平衡にしたときの pH（開放系）
import math
K1, K2, Kw, s, Na = 5.0e-7, 5.0e-11, 1.0e-14, 1.0e-5, 0.10
def bisect(f, lo, hi):          # pH について単調な関数の根
    for _ in range(200):
        mid = (lo + hi)/2
        if f(lo)*f(mid) <= 0: hi = mid
        else: lo = mid
    return (lo + hi)/2
# (1) 閉じた系：近似 [H2CO3]=[CO3^2-]
print('近似 pH', -math.log10(math.sqrt(K1*K2)))
def closed(p):
    h = 10**-p; D = h*h + K1*h + K1*K2
    return Na + h - Na*(K1*h + 2*K1*K2)/D - Kw/h
pc = bisect(closed, 6, 11); h = 10**-pc; D = h*h + K1*h + K1*K2
H2CO3_0 = Na*h*h/D
print('閉鎖系 厳密 pH', pc, ' [H2CO3]', H2CO3_0, ' 大気平衡の', H2CO3_0/s, '倍')
assert abs(pc - 8.30) < 0.02 and 50 < H2CO3_0/s < 150
# (3) 開放系
def opened(p):
    h = 10**-p
    return Na + h - K1*s/h - 2*K1*K2*s/h**2 - Kw/h
po = bisect(opened, 7, 13); h = 10**-po
HCO3, CO3, OH = K1*s/h, K1*K2*s/h**2, Kw/h
print('開放系 pH', po, ' [HCO3-]', HCO3, ' [CO3 2-]', CO3, ' [OH-]', OH)
# 近似：[Na+] = [HCO3-] + 2[CO3^2-]（H+, OH- 無視）→ 2K1K2 s y^2 + K1 s y - Na = 0（y=1/h）
A, B, Cc = 2*K1*K2*s, K1*s, -Na
y = (-B + math.sqrt(B*B - 4*A*Cc))/(2*A); pa = math.log10(y)
print('近似 pH', pa)
assert abs(po - pa) < 0.01 and abs(pa - 10.0) < 1e-6
Cfin = s + HCO3 + CO3; print('炭素の残り', Cfin, ' 失われた割合', (Na - Cfin)/Na)
assert abs((Na - Cfin)/Na - 0.25) < 1e-3
print('OK')
# (5) 0.010 mol/L の場合
Na2 = 0.010
p5 = bisect(lambda p: Na2 + 10**-p - K1*s/10**-p - 2*K1*K2*s/10**(-2*p) - Kw/10**-p, 7, 13)
t = (-1 + math.sqrt(1 + 4*Na2/0.05))/2; print('0.010 mol/L: 厳密 pH', p5, ' 近似 t', t, ' pH', 10 + math.log10(t))
assert abs(p5 - 9.23) < 0.01
print('OK2')
