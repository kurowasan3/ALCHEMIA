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
K1, K2, K3 = 0.40, 1.5, 1.9
x1 = 1/(1+K1); x2 = 1/(1+K2); print(x1, x2); assert sig(x1*100, 2) == 71 and sig(x2*100, 2) == 40
nCO = 1.0 + 1.0/K1; print('CO', nCO); assert abs(nCO-3.5) < 1e-12
y = (K2*2.5 - 1.0)/(1+K2); print('Fe3O4', y); assert abs(y-1.1) < 1e-12 and y > 1/3
def xco(P):
    return (-K3 + math.sqrt(K3**2 + 4*P*K3))/(2*P)
a, b = xco(1.0), xco(3.0); print(a, b); assert sig(a*100, 2) == 72 and sig(b*100, 2) == 54 and a > x1 > b
assert abs(K3**2 + 4*1.0*K3 - 11.21) < 1e-9 and abs(K3**2 + 4*3.0*K3 - 26.41) < 1e-9
dH = 2*(-111) - (-394); print('dH', dH); assert dH == 172
# 各平衡式の元素のつり合い（状態記号を除いて）
for e in ('FeO + CO -> Fe + CO2', 'Fe3O4 + CO -> 3FeO + CO2', 'C + CO2 -> 2CO'):
    assert check_equation(e) is None
print('ALL OK')
