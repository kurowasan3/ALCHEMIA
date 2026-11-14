# THE-ION-038 検証：マレイン酸とフマル酸の中和滴定（近似の破れと第1中和点の明瞭さ）
import json, re, sys
from math import log10, sqrt
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ACIDS = {'maleic': (1.2e-2, 6.0e-7), 'fumaric': (9.0e-4, 4.0e-5)}
def pH(Ka1, Ka2, Vb, Ca=0.100, Va=20.0, Cb=0.100):
    """電荷収支・物質収支を近似なしで解く（二分法）"""
    V = Va + Vb; C = Ca * Va / V; Na = Cb * Vb / V
    def f(p):
        h = 10**-p; D = h * h + Ka1 * h + Ka1 * Ka2
        return h + Na - C * (Ka1 * h + 2 * Ka1 * Ka2) / D - 1e-14 / h
    lo, hi = 0.0, 14.0
    for _ in range(200):
        m = (lo + hi) / 2
        if f(lo) * f(m) <= 0: hi = m
        else: lo = m
    return lo
for k, (a, b) in ACIDS.items():
    print(k, 'pKa', round(-log10(a), 2), round(-log10(b), 2), {v: round(pH(a, b, v), 2) for v in (10, 19, 20, 21, 40)})
assert round(-log10(1.2e-2), 2) == 1.92 and round(-log10(6.0e-7), 2) == 6.22 and round(-log10(9.0e-4), 2) == 3.05 and round(-log10(4.0e-5), 2) == 4.40
# (2) マレイン酸 10 mL：二次方程式（第2段無視）
C = 0.100 * 10 / 30; K1 = 1.2e-2
h = (-(C + K1) + sqrt((C + K1)**2 + 4 * K1 * C)) / 2; print('maleic 10 mL h', h, 'pH', -log10(h))
assert abs(-log10(h) - 2.12) < 0.005 and abs(pH(*ACIDS['maleic'], 10) - 2.12) < 0.01
assert abs(((C + K1)**2 + 4 * K1 * C) * 1e4 - 36.55) < 0.01
# フマル酸 10 mL：pH≈pKa1 がほぼ成り立つ
assert abs(pH(*ACIDS['fumaric'], 10) - 3.05) < 0.05
# (3) 第1中和点（20 mL，C=0.0500）
C1 = 0.0500
for k, (a, b) in ACIDS.items():
    h_simple = sqrt(a * b); h_better = sqrt(a * (b * C1 + 1e-14) / (a + C1)); h_exact = 10**-pH(a, b, 20)
    fr = 1 / (1 + h_exact / a + b / h_exact)
    print(k, 'simple', round(-log10(h_simple), 2), 'better', round(-log10(h_better), 2), 'exact', round(-log10(h_exact), 3), 'HA- frac', round(100 * fr, 1))
a, b = ACIDS['maleic']
assert round(-log10(sqrt(a * b)), 2) == 4.07 and round(-log10(sqrt(a * (b * C1 + 1e-14) / (a + C1))), 2) == 4.12 and round(pH(a, b, 20), 2) == 4.12
assert abs(a * (b * C1) / (a + C1) * 1e9 - 5.81) < 0.01
he = 10**-pH(a, b, 20); assert round(100 / (1 + he / a + b / he), 1) == 98.6
a, b = ACIDS['fumaric']; he = 10**-pH(a, b, 20)
assert round(pH(a, b, 20), 2) == 3.73 and round(100 / (1 + he / a + b / he)) == 70
# (4) 19→21 mL の pH 変化
dM = pH(*ACIDS['maleic'], 21) - pH(*ACIDS['maleic'], 19); dF = pH(*ACIDS['fumaric'], 21) - pH(*ACIDS['fumaric'], 19)
print('jump', round(dM, 2), round(dF, 2)); assert round(dM, 1) == 1.7 and abs(dF - 0.15) < 0.01
assert abs(pH(*ACIDS['fumaric'], 19) - 3.65) < 0.01 and abs(pH(*ACIDS['fumaric'], 21) - 3.80) < 0.01 and abs(pH(*ACIDS['fumaric'], 10) - 3.0) < 0.05
it = json.load(open('/home/claude/chem/out2/items/THE-ION-038/item.json')); k = 0
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') + eq.count('<=>') == 1 and 'boxed' not in eq and not re.search(r'(^|[^A-Z])H?H?A(\^|-| |$)', eq):
            r = check_equation(eq); assert r is None, (eq, r); k += 1; print('eq ok', eq)
print('eq', k)
print('OK')
