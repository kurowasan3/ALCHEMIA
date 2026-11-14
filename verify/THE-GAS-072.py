# THE-GAS-072 検証：クラウジウス・クラペイロンの式による外圧と沸点
import math
R = 8.31
def Tb(P_ratio, T0, dH):
    # log10(P/P0) = -dH/(2.303R)(1/T - 1/T0)
    inv = 1/T0 - 2.303*R*math.log10(P_ratio)/dH
    return 1/inv
lg7 = 0.845 - 1   # log10 0.70 = -0.155
lg2 = 0.301
def Tb_given(lg, T0, dH):
    return 1/(1/T0 - 2.303*R*lg/dH)
for name, T0, dH in (('水', 373, 40.7e3), ('エタノール', 351, 38.6e3)):
    t1 = Tb_given(lg7, T0, dH); t2 = Tb_given(lg2, T0, dH)
    print(name, '0.70e5 Pa:', t1, t1-273, ' 2.0e5 Pa:', t2, t2-273, ' 低下', T0-t1, ' T0^2/dH', T0**2/dH)
    print('   厳密な log', Tb(0.70, T0, dH)-273, Tb(2.0, T0, dH)-273)
# 90℃で沸騰する圧力
lgP = -40.7e3/(2.303*R)*(1/363 - 1/373)
print('log10(P/P0) at 363K', lgP, 'P/P0', 10**lgP, 'P', 1.013e5*10**lgP)
