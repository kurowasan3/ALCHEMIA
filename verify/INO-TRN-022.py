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
Kd, KspI, KspCl, Ka = 1.0e6, 1.0e-12, 2.0e-7, 1.0e-13
# (1) 厳密：c + 2y = 0.0200, y = Kd c^2
c = (-1 + math.sqrt(1 + 8*Kd*0.0200))/(4*Kd); print('Cu+', c); assert sig(c, 2) == 1.0e-4
# (2)
K2 = 1/(Kd*KspI**2); cu2 = 1/(K2*0.10**2); print(K2, cu2); assert sig(K2, 2) == 1.0e18 and sig(cu2, 2) == 1.0e-16
# (3) 厳密に：Cl0 を二分法で。平衡：Cu2+ = 1e-4 のとき [Cl-]^2 = 1/(K' Cu2+)、Cl 収支 Cl0 = [Cl-] + 2*(反応した Cu2+) + 溶けている Cu+ の分は CuCl の Cl とは無関係
Kp = 1/(Kd*KspCl**2); cl = math.sqrt(1/(Kp*1.0e-4)); reacted = 1.0e-2 - 1.0e-4
# 溶けている Cu+ = Ksp/[Cl-]（Cu の粉末から供給されるので Cu2+ 収支には入らないが、CuCl 生成量を正確にする）
cup = KspCl/cl
# Cu 原子収支：Cu2+ 反応量 x と Cu 金属反応量 x → Cu(I) 2x = CuCl + Cu+
nCuCl = 2*reacted - cup
Cl0 = cl + nCuCl; print("K'", Kp, 'Cl eq', cl, 'Cu+', cup, 'Cl0', Cl0); assert sig(Cl0, 2) == 0.040
# 平衡の自己整合（Cu2+ と Cu+ の不均化平衡も成立）
assert abs(1.0e-4/cup**2 - Kd)/Kd < 1e-9
# (4)
K4 = Ka/KspI**2; cu = math.sqrt(5.0e-3/(K4*0.10**4)); print(K4, cu, cu/1e-2); assert sig(K4, 2) == 1.0e11 and sig(cu/1e-2*100, 2) == 0.22
print('ALL OK')
