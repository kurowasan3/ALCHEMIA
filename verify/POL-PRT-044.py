ITEM = '/home/claude/chem5/outh/items/POL-PRT-044/item.json'
# 検証スクリプト（グループ g09）  python3 -I この.py
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert [int(m.group(2)) for m in re.finditer(r'(^|\s)\((\d+)\)', it['answer'])] == list(range(1, np_+1))

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from sympy import symbols, solve, log, Rational as R
c = symbols('c')
bX,KX,bY,KY = R(4),R(2),R('4.5'),R('1.5')
cx = solve(bX-KX*c-1,c)[0]; cy = solve(bY-KY*c-1,c)[0]
print(float(cx), float(cy)); assert round(float(cx),2)==1.50 and round(float(cy),2)==2.33
lo90 = solve(bX-KX*c-0,c)[0]; assert lo90==2
lo99 = solve(bX-KX*c+1,c)[0]; assert lo99==R(5,2) and lo99>cy
import math
def window(d):
    L = math.log10(d)
    lo = (float(bX)+1+L)/float(KX); hi = (float(bY)-1+L)/float(KY)
    return lo, hi
for d in [1,2,5,9.9,10,10.1,20]:
    lo,hi = window(d); print(d, round(lo,4), round(hi,4), lo<=hi+1e-12)
assert window(9.9)[0] > window(9.9)[1] and abs(window(10)[0]-window(10)[1])<1e-12 and window(10)[0]==3.0
# 直接シミュレーション: d=10, c=3.0
d=10; C0=10.0/d; SX=10**(4-2*3.0); SY=10**(4.5-1.5*3.0)
print('X 沈殿率', 1-SX/C0, 'SY', SY, 'C0', C0); assert abs(1-SX/C0-0.99)<1e-9 and SY>=C0-1e-12
# K_X<K_Y なら希釈で窓が狭まる
def win2(d,bX,KX,bY,KY):
    L=math.log10(d); return (bY-1+L)/KY-(bX+1+L)/KX
assert win2(10,4,1.5,4.5,2.0) < win2(1,4,1.5,4.5,2.0)
print('ALL OK')
