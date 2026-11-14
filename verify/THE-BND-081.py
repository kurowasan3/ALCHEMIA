import math
def run(s2, s3, pi):
    out = {}
    for rp in (0.80, 0.50):
        aA = max(2*(rp + 1)/s3, 2.0); aB = max(2*(rp + 1), 2*s2)
        V = 4/3*pi*(rp**3 + 1)
        out[rp] = (aA, aB, (4*V/aB**3)/(V/aA**3))
    return out
for consts in [(1.41, 1.73, 3.14), (math.sqrt(2), math.sqrt(3), math.pi)]:
    o = run(*consts); print(o)
    assert round(o[0.80][0], 2) == 2.08 and o[0.80][0] > 2
    assert round(o[0.80][2], 2) == 0.77 and round(o[0.50][2], 1) == 1.2
    assert round(consts[1] - 1, 2) == 0.73
print('OK')
