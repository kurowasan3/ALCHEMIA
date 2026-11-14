# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-SOL-037/item.json'
from sympy import *
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
print('反応式', len(eqs), '本を確認')

def sig(x, n=2):
    x = float(x)
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def chk(name, got, want, n=2):
    assert sig(got, n) == sig(want, n), (name, float(got), want)
    print('OK', name, sig(got, n))
KD = (2.58e-2 - 3.0e-4)/3.0e-4; chk('(1) KD', KD, 85)
KD = 85
aq_tot = 3.10e-3 - 1.70e-3; i2w = 1.70e-3/KD; i3 = aq_tot - i2w; iminus = 0.100 - i3
K = i3/(i2w*iminus); chk('(2) K', K, 7.0e2); K = 700
Imin3 = (99*KD - 1)/K; chk('(3) [I-]', Imin3, 12); assert Imin3 > 6
Imin4 = (9.9*KD - 1)/K; chk('(4) [I-]eq', Imin4, 1.2)
# 厳密：初濃度 c0 で 99 % 抽出になる c0 を求める
c0, co = symbols('c0 co', positive=True)
ntot = 1.00e-2*0.100
# 有機層 co (0.100 L)，水層 I2 = co/KD，I3- = K*(co/KD)*(c0 - I3-)  (1.00 L)
co_val = 0.01*ntot/0.100   # 有機層に 1 % 残る
w = co_val/KD
i3v = symbols('i3v', positive=True)
i3s = solve(Eq(i3v, K*w*(c0 - i3v)), i3v)[0]
eq = Eq(co_val*0.100 + (w + i3s)*1.00, ntot)
c0v = [s for s in solve(eq, c0) if s > 0][0]
print('厳密な初濃度', float(c0v)); chk('(4) c0 exact', c0v, 1.2)
print('THE-SOL-037 すべてOK')
