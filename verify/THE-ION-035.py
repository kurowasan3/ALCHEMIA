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
from sympy import symbols, nsolve, Rational, log, cbrt
K1, K2, Ksp, Kw = 4.0e-7, 5.0e-11, 5.0e-9, 1.0e-14
K = Ksp*K1/K2; print('K', K); assert abs(K-4.0e-5) < 1e-12
def exact(P):
    """近似なしで電荷収支を解く： 2[Ca]+[H]=[HCO3]+2[CO3]+[OH], [Ca][CO3]=Ksp"""
    c = 3.0e-2*P/1.0e5
    h = symbols('h', positive=True)
    hco3 = K1*c/h; co3 = K2*hco3/h; ca = Ksp/co3
    f = 2*ca + h - hco3 - 2*co3 - Kw/h
    return c, float(nsolve(f, h, 1e-7 if P > 1e3 else 5e-9))
def approx(P):
    c = 3.0e-2*P/1.0e5
    return (K*c/4)**(1/3)
s1 = approx(1.0e5); print('s(1e5)', s1); assert sig(s1, 2) == 6.7e-3
c, h = exact(1.0e5); hco3 = K1*c/h; ca = Ksp/(K2*hco3/h)
print('exact s', ca, 'pH', -float(log(h, 10)))
assert abs(ca - s1)/s1 < 0.002 and round(-float(log(h, 10)), 1) == 6.0
print('CO3/HCO3', K2/h, 'H/HCO3', h/hco3)
# (3)
rem = (1.0e-2 - s1)/1.0e-2; print('remain', rem); assert sig(rem, 2) == 0.33
# (4)
sa, sb = approx(2.0e3), approx(40); dep = (sa - sb)*100; print(sa, sb, dep); assert sig(dep, 2) == 0.13
c, h = exact(2.0e3); sa_ex = Ksp/(K2*(K1*c/h)/h)
c, h = exact(40); sb_ex = Ksp/(K2*(K1*c/h)/h); hco3 = K1*c/h
print('exact', sa_ex, sb_ex, (sa_ex-sb_ex)*100, 'pH cave', -float(log(h, 10)), 'CO3/HCO3', K2/h, 'OH/HCO3', Kw/h/hco3)
assert sig((sa_ex-sb_ex)*100, 2) == 0.13 and round(-float(log(h, 10)), 1) == 8.3
print('ALL OK')
