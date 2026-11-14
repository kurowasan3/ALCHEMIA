# 検証スクリプト（x2 g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x2/out/items/POL-PRT-074/item.json'
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

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 1), (nprob, nans)
print('小問数', nprob, '答え', nans)

from math import log10, sqrt
K1,K2,Kw=5.0e-3,2.0e-10,1.0e-14
assert round(-log10(sqrt(K1*K2)),1)==6.0
def exact(C,Na):
    # 電荷のつり合いを二分法で解く
    lo,hi=1e-14,1.0
    for _ in range(200):
        h=sqrt(lo*hi)
        D=h*h+K1*h+K1*K2
        Ap=C*h*h/D; Am=C*K1*K2/D
        f=Na+h+Ap-Kw/h-Am
        if f>0: hi=h
        else: lo=h
    return -log10(h)
# X: 等電点に合わせたもの（純アラニンの pH はほぼ pI）＋NaOH 0.010
pX=exact(0.100,0.010); pY0=exact(0.100,0.050); pY=exact(0.100,0.060)
print('X後',pX,'Y',pY0,'Y後',pY)
print('純アラニン',exact(0.100,0.0)); assert round(pX,1)==8.7 and round(pY0,1)==9.7 and round(pY,1)==9.9
assert round(9.70-2*0.48,2)==8.74 and round(9.70+0.18,2)==9.88

print('ALL OK')
