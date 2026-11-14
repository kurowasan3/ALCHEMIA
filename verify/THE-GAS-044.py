# THE-GAS-044 検証：可動ピストンの両側に異なる揮発性液体（ベンゼン・エタノール）
# 蒸気圧表（アントワン式から計算した値を 0.1 kPa に丸めたもの）
def pb(t): return 10**(6.90565-1211.033/(t+220.790))*0.133322   # ベンゼン kPa
def pe(t): return 10**(8.20417-1642.89/(t+230.300))*0.133322    # エタノール kPa
table = {t: (round(pb(t),1), round(pe(t),1)) for t in (40, 50, 60, 65, 70, 75, 80)}
print(table)
assert table[40] == (24.4, 17.8) and table[70] == (73.4, 72.2) and table[75] == (86.4, 88.7) and table[60] == (52.2, 46.8)
R = 8.31; V0 = 10.0; nA = 0.120; nB = 0.100
# (1) 40 ℃
T = 313; RT = R*T
pA, pB = table[40]
assert nA > pA*V0/RT and nB > pB*V0/RT          # 最初は両室とも液体が残る
# 両方に液体がある限り pA > pB でつり合わない → A の液体がなくなるまで右へ
# A が全部気体になった後：nA RT / VA = pB（B に液体が残る）
VA = nA*RT/pB; VB = 2*V0 - VA
print('(1) VA =', VA, 'VB =', VB)
assert VA < 2*V0
nBv = pB*VB/RT
assert nBv < nB and pB < pA                       # B に液体が残り，A は飽和に達しない
print('    B の蒸気', nBv, '液体', nB - nBv)
assert round(VA,1) == 17.5 and round(nB-nBv,4) == 0.0831
# (2) 蒸気圧が等しい温度（表の直線補間）
t1, t2 = 70, 75
d1 = table[t1][0]-table[t1][1]; d2 = table[t2][0]-table[t2][1]
tx = t1 + (t2-t1)*d1/(d1-d2)
print('(2) 交点', tx, ' アントワン式で直接', [t/100 for t in range(7000, 7500) if pb(t/100) <= pe(t/100)][0])
assert 71 < tx < 72.5
# (3) 60 ℃
T = 333; RT = R*T
p = (nA+nB)*RT/(2*V0); VA3 = 2*V0*nA/(nA+nB)
print('(3) p =', p, 'VA =', VA3)
assert p < min(table[60])
assert round(VA3,1) == 10.9
print('OK')
