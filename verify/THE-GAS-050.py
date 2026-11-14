"""THE-GAS-050 検証"""
P = lambda h: 1.00e5 + 1.00e4 * h
rho = P(20.0) * 44.0 / (8.31e3 * 283)
ratio = P(20.0) / P(0)
h = (4.50e6 - 1.00e5) / 1.00e4
print(rho, ratio, h)
assert f'{rho:.2g}' == '5.6' and abs(ratio - 3.0) < 1e-12 and f'{h:.2g}' == '4.4e+02'
assert 35 > 31   # 臨界温度より高い
print('THE-GAS-050 OK')
