# THE-EQL-163 検証：大気中の硫黄化合物の寿命（擬一次反応・競合する酸化剤）
import math
P, R, T, NA = 1.013e5, 8.31, 298, 6.02e23
ntot = P/(R*T)*NA/1e6          # 個/cm3
nO3 = 30e-9*ntot
print('air', ntot, 'O3', nO3)
assert f'{nO3:.3g}' == '7.39e+11'
OH = 1.0e6
k = {'DMS': (4.8e-12, 1.0e-18), 'H2S': (4.7e-12, 0.0), 'COS': (2.0e-15, 0.0)}
tau = {}
for X, (kOH, kO3) in k.items():
    a, b = kOH*OH, kO3*nO3
    tau[X] = 1/(a + b)
    print(X, 'kOH[OH]', a, 'kO3[O3]', b, 'frac OH', a/(a+b), 'tau s', tau[X], 'days', tau[X]/86400, 'years', tau[X]/86400/365)
assert f"{tau['DMS']/86400:.2g}" == '2.1'
assert f"{tau['H2S']/86400:.2g}" == '2.5'
assert f"{tau['COS']/86400/365:.2g}" == '16'
fr = 4.8e-6/(4.8e-6 + 1.0e-18*nO3)
assert f'{fr*100:.2g}' == '87'
surv = math.exp(-365*86400/tau['COS'])
print('COS survive 1 yr', surv, 'H2S', math.exp(-365*86400/tau['H2S']))
assert f'{surv:.2g}' == '0.94'
print('half-life COS yr', tau['COS']*math.log(2)/86400/365)
# 100 ppt DMS と OH・O3 の比
nDMS = 100e-12*ntot
print('DMS', nDMS, 'DMS/OH', nDMS/OH, 'DMS/O3', nDMS/nO3)
