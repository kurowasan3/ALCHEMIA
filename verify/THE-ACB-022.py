# 検証スクリプト（共通ヘッダ＋問題固有の検証）  python3 -I <this file>
import json, os, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
ID = os.path.splitext(os.path.basename(__file__))[0]
ITEM = json.load(open(f'/home/claude/chem/out2/items/{ID}/item.json', encoding='utf-8'))
def _ce(s):
    out = []; i = 0
    while (i := s.find('\\ce{', i)) >= 0:
        j = i + 4; dep = 1
        while dep:
            dep += {'{': 1, '}': -1}.get(s[j], 0); j += 1
        out.append(s[i+4:j-1]); i = j
    return out
EQS = []
for _k in ('problem', 'answer', 'solution'):
    for _c in _ce(ITEM[_k]):
        if len(re.findall(r'<=>|->|<-', _c)) == 1 and not re.search(r'\\boxed|\\phantom|\\text|\bR\b|\bn\b', _c):
            EQS.append(_c)
for _c in EQS:
    _r = check_equation(_c)
    assert _r is None, (_c, _r)
print(f'{ID}: 反応式 {len(EQS)} 本はすべて原子数・電荷がつり合っている')
# 小問数と答えの番号の一致
_nq = ITEM['problem'].count('\\item')
_na = len(re.findall(r'(?:^|\s)\((\d+)\)', ITEM['answer']))
assert _nq == _na or (_nq == 0 and _na == 0), (_nq, _na)
print(f'{ID}: 小問 {_nq} 個と答え {_na} 個が対応')
M = dict(H=1.0, C=12.0, N=14.0, O=16.0, Na=23.0, Mg=24.0, Al=27.0, Si=28.0, P=31.0, S=32.0, Cl=35.5, K=39.0,
         Ca=40.0, Ti=48.0, Cr=52.0, Mn=55.0, Fe=56.0, Co=59.0, Cu=63.5, Zn=65.4, Br=80.0, Sr=88.0, Ag=108.0, I=127.0, Ba=137.0, Pt=195.0)
def mass(f, masses=None):
    m = masses or M
    c, _ = parse_formula(f)
    return sum(m[e] * float(n) for e, n in c.items())
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
# ---- 問題固有の検証 ----
from sympy import symbols, nsolve, sqrt, log, solve
K1, K2, K3, Kw = 7.0e-3, 6.0e-8, 4.0e-13, 1.0e-14
nP = 0.100*3.00e-3; nH = 0.100*(8.00-3.00)*1e-3
assert abs(nP/0.0100 - 0.0300) < 1e-12 and abs(nH/0.0100 - 0.0500) < 1e-12
# 完全な電荷収支で滴定曲線を計算（第1・第2当量点の pH，第3当量点付近の傾き）
def pH_at(vNaOH_mL):
    """電荷収支を pH について二分法で解く（数値的に厳密）"""
    import math
    V = (10.0 + vNaOH_mL)/1000
    C = nP/V; Cl = nH/V; Na = 0.100*vNaOH_mL/1000/V
    def f(pH):
        h = 10**(-pH)
        D = h**3 + K1*h**2 + K1*K2*h + K1*K2*K3
        a1 = K1*h**2/D; a2 = K1*K2*h/D; a3 = K1*K2*K3/D
        return Na + h - Cl - Kw/h - C*(a1 + 2*a2 + 3*a3)
    lo, hi = 0.0, 14.0
    for _ in range(200):
        mid = (lo+hi)/2
        if f(lo)*f(mid) <= 0: hi = mid
        else: lo = mid
    return (lo+hi)/2
for v in (7.9, 8.0, 8.1, 10.9, 11.0, 11.1, 13.5, 14.0, 14.5):
    print(v, round(pH_at(v), 3))
# 第2当量点：完全解と (c) の式が一致すること
p2 = pH_at(11.0); print('pH(eq2) exact', p2)
C = nP/0.0210
hc = ((K2*(K3*C + Kw))/C)**0.5; pc = -float(log(hc, 10)); print('pH (c)', pc)
assert round(pc, 1) == 9.6 and abs(p2 - pc) < 0.03
pa = -float(log((K2*K3)**0.5, 10)); print('pH (a)', pa); assert round(pa, 2) == 9.81
h = (K2*K3)**0.5; print('PO4', K3*C/h, 'OH', Kw/h); assert sig(K3*C/h, 2) == 3.7e-5 and sig(Kw/h, 2) == 6.5e-5
# (2) 加水分解の割合
a = symbols('a', positive=True); Cp = nP/0.0240; Kh = Kw/K3
al = [float(r) for r in solve(Cp*a**2 - Kh*(1-a), a) if 0 < float(r) < 1][0]; print('alpha', al); assert sig(al, 2) == 0.73
# 第3当量点付近の pH 変化が小さいこと（13.5→14.5 mL）と第2当量点付近（10.9→11.1）の比較
d3 = pH_at(14.5) - pH_at(13.5); d2 = pH_at(11.1) - pH_at(10.9); print('jump eq2 (0.2mL)', d2, 'eq3 (1.0mL)', d3)
assert d3 < 0.5 and d2 > d3
print('ALL OK')
