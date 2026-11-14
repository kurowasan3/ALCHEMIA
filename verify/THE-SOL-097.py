# THE-SOL-097 検証：不純物 O2 を含む窒素での繰り返しの脱気（ヘンリーの法則・漸化式の不動点）
from fractions import Fraction as F
H = F(30)        # mL（25℃,1.00e5 Pa 換算）/ 水 1.00 L / O2 分圧 1.00e5 Pa
V0 = H * F(1, 5) # 空気（O2 20.0%）と平衡：6.00 mL
def run(eps, vg, nmax=50):
    seq = [V0]
    for _ in range(nmax):
        tot = seq[-1] + vg * eps           # 水中の O2 + 注入気体中の O2（mL）
        p = tot / (vg + H)                 # O2 の分圧（1.00e5 Pa 単位）
        seq.append(H * p)
    return seq
eps = F(1, 10**4); vg = F(150)
s = run(eps, vg)
print('V1', float(s[1]))
assert f'{float(s[1]):.3g}' == '1'
Vstar = H * eps
print('fixed point', float(Vstar))
assert Vstar == F(3, 1000)
for a, b in zip(s, s[1:]):
    assert b - Vstar == (a - Vstar) * H / (vg + H)
n = next(k for k, v in enumerate(s) if v <= V0 / 1000)
s0 = run(F(0), vg)
n0 = next(k for k, v in enumerate(s0) if v <= V0 / 1000)
print('n with impurity', n, 'pure', n0)
assert (n, n0) == (5, 4)
assert all(v > V0 / 10000 for v in s)       # 1/10000 には届かない
# 1/10000 に（極限で）届く条件：30 eps < 6e-4 → eps < 2e-5
assert H * F(2, 10**5) == V0 / 10000
# 気相 0.300 L の場合
s3 = run(eps, F(300))
n3 = next(k for k, v in enumerate(s3) if v <= V0 / 1000)
print('n (0.300 L)', n3, 'limit', float(s3[-1]))
assert n3 == 4 and abs(s3[-1] - Vstar) < F(1, 10**12)
print('OK')
