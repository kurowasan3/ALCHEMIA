# 検証スクリプト（x2 g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x2/out/items/ORG-HYD-054/item.json'
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

H,C,Br,Ag=1.0,12.0,79.9,107.9
mB=3*C+3*H+Ag; mD=3*C+6*H+2*Br; mAg2=2*C+2*Ag; mBr2=2*Br
ny=0.0600; ne=0.0800; na=0.0600
assert round(ny*mB,2)==8.81 and round((0.200-ny)*22.4,2)==3.14 and round(ne*mD,1)==16.1
print('比', na/0.02, ne/0.02, ny/0.02)
assert abs(round(16.1/mD,4)-0.0798)<1e-4
b=(ne+2*ny)*mBr2; print('Br2',b); assert round(b,1)==32.0
s=ny*mAg2; print('Ag2C2',s); assert round(s,1)==14.4

print('ALL OK')
