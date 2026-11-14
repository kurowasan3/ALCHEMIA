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
import math
K, KBa, KSr = 3.0e14, 1.2e-10, 2.2e-5
Kp = 1/(K*KBa**2); print("K'", Kp); assert sig(Kp, 2) == 2.3e5
# (2)
h = math.sqrt(Kp*(1e-2)**2*1e-2); pH = -math.log10(h); print('h', h, 'pH', pH); assert round(pH, 1) == 0.3
# 厳密：沈殿しはじめの条件を全 Cr 収支込みで
c = KBa/1e-2
def f(hh): return c + 2*K*c**2*hh**2 - 2.0e-2
lo, hi = 0.01, 10
for _ in range(200):
    m = (lo*hi)**0.5
    if f(lo)*f(m) <= 0: hi = m
    else: lo = m
print('exact threshold pH', -math.log10(lo)); assert round(-math.log10(lo), 1) == 0.3
assert -math.log10(10**-4.7) > 0.3
# (3) 厳密解（二分法）：p = 反応した Cr2O7
def g(p):
    Ba = 3.0e-2 - 2*p; H = 2*p
    cr = K*(KBa/Ba)**2*H**2         # 平衡の Cr2O7
    return (1.0e-2 - p) - cr - (KBa/Ba)/2
lo, hi = 1e-4, 0.999e-2
for _ in range(200):
    m = (lo+hi)/2
    if g(lo)*g(m) <= 0: hi = m
    else: lo = m
p = lo; rem = (1.0e-2 - p)/1.0e-2; print('remain exact', rem)
H = 2.0e-2; Ba = 1.0e-2; cr = H**2/(Kp*Ba**2); rem2 = 2*cr/2.0e-2; print('approx', cr, rem2)
assert sig(rem*100, 2) == 0.17 and sig(rem2*100, 2) == 0.17
# (4)
a = 2*K*(1e-5)**2; cc = (-1 + math.sqrt(1 + 4*a*1.9e-2))/(2*a); print('CrO4', cc, 1+4*a*1.9e-2); assert sig(cc, 2) == 5.5e-4
ba_rem = KBa/cc/1e-3; print('Ba remain', ba_rem); assert sig(ba_rem*100, 2) == 0.022
assert 1e-2*cc < KSr
a7 = 2*K*(1e-7)**2; c7 = (-1 + math.sqrt(1 + 4*a7*1.9e-2))/(2*a7); print('pH7 CrO4', c7, 1e-2*c7); assert 1e-2*c7 > KSr and sig(c7, 2) == 1.7e-2
# pH7 で Sr が沈殿するとき Ba の分は無視できない？ → 判断（沈殿するか否か）には影響しない
print('ALL OK')
