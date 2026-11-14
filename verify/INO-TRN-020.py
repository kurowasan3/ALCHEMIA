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
assert len(EQS) >= 8, EQS
def ox_Mn(f):
    c, q = parse_formula(f)
    return (F(q) - sum(F({'O': -2, 'H': 1}[e]) * n for e, n in c.items() if e != 'Mn')) / c['Mn']
v = {f: ox_Mn(f) for f in ('MnO4-', 'MnO4^2-', 'MnO2', 'Mn^2+')}
print(v); assert list(v.values()) == [7, 6, 4, 2]
e = {k: 7 - v[f] for k, f in (('X', 'Mn^2+'), ('Y', 'MnO2'), ('Z', 'MnO4^2-'))}
so3 = {k: F(1, 100) * n / 2 for k, n in e.items()}
print(so3); assert so3 == {'X': F(25, 1000), 'Y': F(15, 1000), 'Z': F(5, 1000)}
# 係数比とも一致
for c, k in (('2MnO4- + 5SO3^2-', 'X'), ('2MnO4- + 3SO3^2-', 'Y'), ('2MnO4- + SO3^2-', 'Z')):
    assert any(eq.startswith(c) for eq in EQS), c
print('ALL OK')
