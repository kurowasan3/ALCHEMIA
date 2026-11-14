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
from fractions import Fraction as F
nFe2 = 0.0200*17.0e-3*5; nFe = 0.0200*19.0e-3*5; nFe3 = nFe - nFe2
print(nFe2, nFe, nFe3); assert abs(nFe2-1.70e-3) < 1e-12 and abs(nFe3-2.0e-4) < 1e-12
mO = 0.1384 - nFe*56; nO = mO/16; z1 = nFe/nO
nO2 = (2*nFe2 + 3*nFe3)/2; z2 = nFe/nO2
print('z (mass)', z1, 'z (charge)', z2); assert round(z1, 3) == 0.95 and round(z2, 3) == 0.95
print('Fe3+ frac', nFe3/nFe); assert sig(nFe3/nFe, 2) == 0.11
NA = 6.0e23; V = 0.430**3*1e-21   # nm^3 -> cm^3
V_given = 0.0795e-21
dA = (4*0.95*56 + 4*16)/(NA*V_given); dB = (4*56 + 4/0.95*16)/(NA*V_given)
dA2 = (4*0.95*56 + 4*16)/(NA*V); dB2 = (4*56 + 4/0.95*16)/(NA*V)
print(dA, dB, dA2, dB2); assert sig(dA, 3) == 5.80 and sig(dB, 3) == 6.11 and sig(dA2, 3) == 5.80 and sig(dB2, 3) == 6.11
vac = 4*(1-0.95); fe3 = 4*0.95*(0.10/0.95); print(vac, fe3); assert abs(vac-0.20) < 1e-9 and abs(fe3-0.40) < 1e-9
# 酸化数の確認（Fe_zO の平均酸化数）
assert abs((2*(0.95-0.10) + 3*0.10) - 2) < 1e-12
print('ALL OK')
