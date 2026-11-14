"""ORG-DET-030 検証：2 つのアミンの酸の奪い合いと分配平衡（厳密解と近似解の照合）"""
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from mpmath import mp, mpf, sqrt, log10
mp.dps = 30
KaA, KaB, KDA, KDB = mpf('2.5e-5'), mpf('5.0e-10'), 10, 2
V = mpf('0.1')
def solve(nA, nB, nH, Vw=V, Ve=V):
    def parts(h):
        Aw = nA/(Vw*(1+h/KaA)+KDA*Ve); Bw = nB/(Vw*(1+h/KaB)+KDB*Ve)
        return dict(AH=Aw*h/KaA*Vw, BH=Bw*h/KaB*Vw, Aw=Aw*Vw, Bw=Bw*Vw)
    f = lambda h: (lambda p: p['AH']+p['BH']+h*Vw-nH)(parts(h))
    lo, hi = mpf('1e-14'), mpf('1')
    for _ in range(300):
        mid = sqrt(lo*hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    p = parts(lo); p['h'] = lo; return p
K = KaB/KaA; print('K', K); assert abs(K-mpf('2e-5')) < 1e-12
n = mpf('0.01')
r = solve(n, n, n)
# 近似解
q = sqrt(3*K/11); x = n*q/(1+q); print('x exact', r['AH'], 'approx', x)
assert abs(r['AH']-x)/x < 0.01 and round(float(x), 6) == 2.3e-5
assert r['h']*V < x/100                      # 水層の H+ は無視できる
totA = r['AH']+r['Aw']; print('aq aniline', totA, 100*totA/n)
assert round(float(totA), 5) == 9.3e-4 and round(float(100*totA/n), 1) == 9.3
assert r['AH']/totA < 0.03
print('B free aq', r['Bw']); assert abs(r['Bw']-mpf('7.8e-6')) < mpf('1e-7')
print('[H+]', r['h']); assert abs(r['h']-mpf('6.4e-7')) < mpf('0.1e-7')
# 洗浄
r2 = solve(totA, r['BH']+r['Bw'], n)
exact_ratio = (r2['AH']+r2['Aw'])/totA
approx_ratio = (r['AH']+r['Aw']/11)/totA
print('wash ratio exact', exact_ratio, 'approx', approx_ratio)
assert round(float(exact_ratio), 1) == 0.1 and round(float(approx_ratio), 1) == 0.1
# 10% 過剰
r3 = solve(n, n, mpf('0.011'))
t3 = r3['AH']+r3['Aw']; print('excess aq aniline', t3, 100*t3/n, 'pH', -log10(r3['h']))
assert round(float(100*t3/n)) == 18
approx3 = mpf('1.0e-3') + (n-mpf('1.0e-3'))/11; assert abs(approx3-t3)/t3 < 0.01
# 不足ぎみ
r4 = solve(n, n, mpf('0.009')); print('deficit aq aniline %', 100*(r4['AH']+r4['Aw'])/n, 'B recovered %', 100*(r4['BH']+r4['Bw'])/n)
it = json.load(open('/home/claude/chem5/outh/items/ORG-DET-030/item.json'))
for key in ('problem', 'answer', 'solution'):
    for m in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[key]):
        if re.search(r'->|<=>', m):
            res = check_equation(m); print(m, res or 'OK'); assert res is None
print('ALL OK')
