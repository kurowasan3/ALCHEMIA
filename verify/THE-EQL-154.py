# 検証スクリプト THE-EQL-154   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-EQL-154/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import sympy as sp
kf, K, ki, kj, O3, NO2 = 3.5e-17, 3.0e-11, 2.0e-4, 5.0e-3, 1.0e12, 2.5e11
# 定常状態を記号で解いて式を確認
x, y = sp.symbols('NO3 N2O5', positive=True)
Kf, KK, Ki, Kj, o, n = sp.symbols('k_f K k_i k_j O3 NO2', positive=True)
sol = sp.solve([sp.Eq(y, KK*n*x), sp.Eq(Kf*n*o, Ki*y + Kj*x)], [x, y], dict=True)[0]
assert sp.simplify(sol[x] - Kf*n*o/(Kj + Ki*KK*n)) == 0
# 完全な速度式（(g),(h) の速度定数 kg, kh を有限にして平衡が速い極限）で数値的に確認
kh = 50.0; kg = K*kh   # K = kg/kh
from scipy.optimize import fsolve
def eqs(z):
    a, b = z  # NO3, N2O5（NO2,O3 は一定）
    return [kf*NO2*O3 - kg*NO2*a + kh*b - kj*a, kg*NO2*a - kh*b - ki*b]
a, b = fsolve(eqs, [1e9, 1e10])
NO3 = kf*NO2*O3/(kj + ki*K*NO2); N2O5 = K*NO2*NO3
print('定常', NO3, N2O5, '完全式', a, b)
assert abs(a/NO3 - 1) < 1e-2 and abs(b/N2O5 - 1) < 1e-2
vf, vi, vj = kf*NO2*O3, ki*N2O5, kj*NO3
vN = 2*vi + vj; print('(3)', NO3, N2O5, vN, vN/vf)
assert (sig(NO3,2), sig(N2O5,2), sig(vN,2), sig(vN/vf,2)) == (1.3e9, 1.0e10, 1.1e7, 1.2)
assert abs(vi + vj - vf) < 1e-6*vf
N = NO2 + NO3 + 2*N2O5; tau = N/vN/3600; print('(4)', N, tau); assert sig(tau,2) == 7.0
# (5) NO2 を大きくしたときの極限
for big in [1e13, 1e14, 1e15]:
    n3 = kf*big*O3/(kj + ki*K*big); n5 = K*big*n3
    t = (big + n3 + 2*n5)/(2*ki*n5 + kj*n3)
    print(big, t)
tl = 1/(2*kf*O3) + 1/ki; print('(5)', tl, tl/3600); assert sig(tl/3600,2) == 5.4 and abs(t/tl - 1) < 1e-2
# 平衡定数の妥当性（298 K で K ≈ 2.9e-11 cm^3 程度）
for e in ['NO2 + O3 -> NO3 + O2', 'NO2 + NO3 -> N2O5', 'N2O5 -> NO2 + NO3', 'N2O5 + H2O -> 2HNO3']:
    assert check_equation(e) is None

print('検証OK')
