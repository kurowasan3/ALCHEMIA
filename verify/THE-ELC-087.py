ITEM = '/home/claude/chemk/out/items/THE-ELC-087/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import numpy as np
lc = np.array([-4,-3,-2,-1.0]); EK = np.array([-34.1,22.9,79.9,136.9]); ECl = np.array([29.0,-27.0,-83.0,-139.0])
assert np.allclose(np.diff(EK), 57.0) and np.allclose(np.diff(ECl), -56.0)
sK, iK = np.polyfit(lc, EK, 1); sC, iC = np.polyfit(lc, ECl, 1); print(sK, iK, sC, iC)
assert abs(iK-193.9) < 1e-6 and abs(iC+195.0) < 1e-6
K = 10**((27.5-iK)/sK); Cl = 10**((-65.1-iC)/sC); Na = Cl-K
print(K, Cl, Na)
assert sig(K,3) == 1.20e-3 and sig(Cl,3) == 4.79e-3 and sig(Na,3) == 3.59e-3
mK = K*(39.1+35.5); mN = Na*(23.0+35.5); print(mK, mN, mK+mN, mK/0.300*100)
assert sig(mK,3) == 0.0898 and sig(mN,3) == 0.210 and sig(mK+mN,3) == 0.300 and sig(mK/0.300*100,2) == 30
# 設定値からの逆算
k0 = 0.0900/74.6; n0 = 0.210/58.5
print(iK+sK*np.log10(k0), iC+sC*np.log10(k0+n0))
assert abs(iK+sK*np.log10(k0)-27.5) < 0.06 and abs(iC+sC*np.log10(k0+n0)+65.1) < 0.06

print('全検証OK:', ITEM)
