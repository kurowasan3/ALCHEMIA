ITEM = '/home/claude/chemx/tokodai/outh/items/THE-BND-083/item.json'
# 検証スクリプト（東工大 g03）  python3 -I この.py
# 共通の検証ヘッダ（各検証スクリプトにそのまま埋め込む）
import sys, json, re, math
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
nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')
import numpy as np
rA,rB,rH=1.25,1.45,0.40; s2,s3=1.41,1.73
a0=2*(rA+rB)/s3; chk('a0',a0,3.12,3)
# 幾何：B 頂点・A 中心。P=(a/2,a/2,0) Q=(a/2,0,0)
def dists(site,a):
    import itertools
    B=[np.array(p)*a for p in itertools.product([0,1],repeat=3)]
    A=[np.array([0.5,0.5,0.5])*a+np.array(t)*a for t in itertools.product([-1,0],repeat=3)]
    s=np.array(site)*a
    dA=sorted(np.linalg.norm(x-s) for x in A); dB=sorted(np.linalg.norm(x-s) for x in B)
    return dA,dB
dA,dB=dists([0.5,0.5,0],1.0); print('P', dA[:2], dB[:4]); assert abs(dA[0]-0.5)<1e-9 and abs(dA[1]-0.5)<1e-9 and abs(dB[3]-2**-0.5)<1e-9 and dA[2]>0.8
dA,dB=dists([0.5,0,0],1.0); print('Q', dA[:4], dB[:2]); assert abs(dB[1]-0.5)<1e-9 and abs(dA[3]-2**-0.5)<1e-9 and dB[2]>1.0
assert a0/2 < rA+rH and a0/s2 > rB+rH
aP=max(2*(rA+rH), s2*(rB+rH), a0); aQ=max(2*(rB+rH), s2*(rA+rH), a0)
chk('aP',aP,3.30,3); chk('aQ',aQ,3.70,3)
# 膨張後に A-B が重ならない
assert aP*s3/2 >= rA+rB
chk('vol',(3.30/3.12)**3,1.2,2); chk('vol2',(aP/a0)**3,1.2,2)
chk('H%',3/(59+48+3)*100,2.7,2)
dens=3*1.0/(6.0e23*(3.30e-8)**3); chk('dens',dens,0.14,2); chk('ratio',dens/0.071,2.0,2)
print('ALL OK')
