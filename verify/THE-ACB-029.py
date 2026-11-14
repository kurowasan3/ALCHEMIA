ITEM = '/home/claude/chem4/out/items/THE-ACB-029/item.json'
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
R=sp.Rational
c=R(200,10000)*12/2/10; assert c==R(12,1000)
n0=c*R(1,10); left=R(200,10000)*R(26,1000)/2*4; co2=n0-left; print(co2)
assert co2==R(16,100000)
pct=float(co2*R(224,10)/10*100); print(pct); assert sig(pct,2)==0.036
# (4) 沈殿ごと MO：滴定値は全 Ba 量 → CO2 計算値 0
# (5) NaOH
oh0=2*n0; oh=oh0-2*co2; co3=co2
hcl25=(oh+co3)/4; vol=float(hcl25/R(200,10000)*1000); print('NaOH titr mL',vol); assert abs(vol-28.0)<1e-9
co2_app=(oh0-4*hcl25)/2; print(co2_app, co2_app/co2); assert co2_app/co2==R(1,2)
print('全検証OK:', ITEM)
