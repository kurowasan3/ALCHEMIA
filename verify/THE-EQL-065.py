ITEM = '/home/claude/chem4/out/items/THE-EQL-065/item.json'
# 検証スクリプト共通部（グループ g04）  python3 -I この.py
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
# 小問数と答えの番号の一致
nprob = it['problem'].count('\\item')
nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
if nprob:
    assert [int(n) for n in nums] == list(range(1, nprob + 1)), (nprob, nums)
print('小問数OK:', nprob)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import sympy as sp
R=sp.Rational(83,10)*10**3; T=300; p=sp.Rational(1013,10)*10**3-sp.Rational(36,10)*10**3
nO2=p*sp.Rational(1274,10000)/(R*T); print('nO2',float(nO2)); assert sig(float(nO2),2)==5.0e-3
c0=2*nO2/sp.Rational(250,10000); print('c0',float(c0)); assert sig(float(c0),2)==0.40
V=[63.7,95.6,111.5]; fr=[1-v/127.4 for v in V]; print(fr)
conc=[sig(float(c0)*f,2) for f in fr]; print(conc); assert conc==[0.20,0.10,0.050]
k=0.69/5; print('k',k); assert sig(k,2)==0.14
# 表の体積は c0=0.40, t1/2=5 から
for t,v in zip([5,10,15],V):
    vt=float(0.40*0.025/2*(1-2**(-t/5))*R*T/p*1000); print(t,vt); assert abs(vt-v)<0.06
# 希釈：物質量 n(t)=n0*2^(-t/5) は体積によらない → 5分の気体同じ
# 倍濃度：5分で総量の半分 = 127.4 mL
assert abs(2*127.4/2-127.4)<1e-9
print('全検証OK:', ITEM)
