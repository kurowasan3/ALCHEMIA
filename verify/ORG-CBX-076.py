# ORG-CBX-076 verification: palmitic acid monolayer
NA = 6.0e23
n = 2.56e-3 / 256 * 0.50 / 100; N = n * NA; A = 60e-4 / N
print('n', n, 'N', N, 'A', A)
assert abs(n - 5.0e-8) < 1e-12 and abs(N - 3.0e16) < 1e12 and abs(A - 2.0e-19) < 1e-22
v = 256 / (0.85 * NA) * 1e-6; t = v / A; s = A ** 0.5
print('v', v, 't', t, 'ratio1', t / s)
assert round(t * 1e9, 1) == 2.5 and round(t / s, 1) == 5.6
t2 = t / 2; s2 = (2 * A) ** 0.5; print('ratio2', t2 / s2); assert round(t2 / s2, 1) == 2.0
assert abs(20 ** 0.5 - 4.47) < 0.005 and abs(40 ** 0.5 - 6.32) < 0.005
# molar masses
M = lambda C, H, O: 12 * C + 1.0 * H + 16 * O
assert M(16, 32, 2) == 256 and M(22, 44, 2) == 340
tb = 340 / (0.85 * NA) * 1e-6 / A; print('behenic thickness', tb); assert round(tb * 1e9, 1) == 3.3
print('OK')
