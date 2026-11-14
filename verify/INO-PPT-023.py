# 検証スクリプト（自動生成ヘッダ＋問題固有の検証）  python3 -I <this file>
import json, os, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
ID = os.path.splitext(os.path.basename(__file__))[0]
ITEM = json.load(open(f'/home/claude/chem/out/items/{ID}/item.json', encoding='utf-8'))
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
M = dict(H=1.0, He=4.0, C=12.0, N=14.0, O=16.0, Na=23.0, Mg=24.0, Al=27.0, Si=28.0, P=31.0, S=32.0, Cl=35.5, K=39.0,
         Ca=40.0, Ti=48.0, Cr=52.0, Mn=55.0, Fe=56.0, Co=59.0, Cu=63.5, Zn=65.4, Ga=69.7, Br=80.0, Ag=108.0, I=127.0, Pt=195.0, Au=197.0, Pb=207.0)
def mass(f, masses=None):
    m = masses or M
    c, _ = parse_formula(f)
    return sum(m[e] * float(n) for e, n in c.items())
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
# ---- 問題固有の検証 ----
from fractions import Fraction as F
assert len(EQS) >= 4, EQS
assert any(e.startswith('2CoCl2 + 2NH4Cl + 10NH3 + H2O2') for e in EQS)
# Co の酸化数：[Co(NH3)6]^3+ では NH3 は中性
c, q = parse_formula('[Co(NH3)6]^3+'); assert q == 3
c2, q2 = parse_formula('CoCl2'); assert q2 == 0 and c2['Cl'] == 2   # Co = +2
h2o2 = F(2, 100) / 2; print(h2o2); assert h2o2 == F(1, 100)
print('ALL OK')
