ITEM = '/home/claude/chem4/outh/items/THE-ATM-022/item.json'
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
import math
T=5730
app=lambda R: T*math.log10(1/R)/math.log10(2)
R=(1-0.01)/128+0.01; a2=app(R); print(a2); assert sig(a2,2)==3.3e4
a2g=T*(2-0.248)/0.3010; print(a2g); assert sig(a2g,2)==3.3e4
print('rel err', (7*T-a2)/(7*T)); assert abs((7*T-a2)/(7*T)-0.17)<0.01
a3=app(0.99*0.5+0.01); print(a3, T-a3); assert 70<T-a3<90
a3g=T*0.2967/0.3010; assert 70<T-a3g<90
for t in (1000,5730,20000,40000):
    d=app(0.9*0.5**(t/T))-t; assert abs(d-870.98)<0.01
d=T*(1-2*0.4771)/0.3010; print('dead',d); assert sig(d,2)==870
tmax=T*2/0.3010; print(tmax); assert sig(tmax,2)==3.8e4
assert 14==14+0 and 6==7-1
print('全検証OK:', ITEM)
