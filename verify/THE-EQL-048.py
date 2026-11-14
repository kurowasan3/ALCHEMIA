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
Ka, Kb = 1.0e-8, 1.0e20
x = math.sqrt(Ka*1.0**2); conv = 2*x/(1.0 + 2*x); print(x, conv); assert sig(x, 2) == 1.0e-4 and sig(conv*100, 2) == 0.020
K = Ka*Kb; assert sig(K, 2) == 1.0e12
# 厳密：全圧 1.0、Cl2 が y だけ残る。TiCl4 = t, CO = 2t, O2 = q, y + t + 2t + q = 1
# Cl 収支: 2*y_in = 2*y + 4t -> 単位量の Cl2 を流したとして比で解く：p_TiCl4 p_CO^2 / p_Cl2^2 = K
def f(t):
    y = 1 - 3*t
    return t*(2*t)**2 - K*y**2
lo, hi = 0.3, 1/3
for _ in range(200):
    m = (lo+hi)/2
    if f(lo)*f(m) <= 0: hi = m
    else: lo = m
t = lo; y = 1 - 3*t; print('pCl2 exact', y); assert sig(y, 2) == 3.8e-7
print('pO2', (2*t)**2/Kb); assert (2*t)**2/Kb < 1e-20
assert abs((1/3)*(2/3)**2/1e12 - 1.48e-13) < 1e-15
dHa = -763 - (-944); dH = (-763 + 2*(-111)) - (-944); print(dHa, dH); assert dHa == 181 and dH == -41
for e in ('TiO2 + 2Cl2 -> TiCl4 + O2', '2C + O2 -> 2CO', 'TiO2 + 2C + 2Cl2 -> TiCl4 + 2CO'):
    assert check_equation(e) is None
print('ALL OK')
