ITEM = '/home/claude/chem6/outh/items/ORG-STR-050/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if should_check(b):
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
print('検査した式の数:', n_eq)

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


from fractions import Fraction as F
valid=[(2,4),(2,2),(3,6),(3,4),(4,8),(4,6),(4,4),(4,2)]
out=[]
for n in range(1,5):
    for m,h in valid:
        k=(2*m+2-h)//2
        VQ=F(15,k); VP=10-VQ
        if VP<=0: continue
        if n*VP+m*VQ==35 and (n+1)*VP+F(h,2)*VQ==30: out.append((n,m,h,k,VP,VQ))
print(out)
assert [(o[0],o[1],o[2],o[4],o[5]) for o in out]==[(2,4,6,F(5,2),F(15,2)),(3,4,4,5,5)]
for n,m,h,k,VP,VQ in out:
    C=n*VP+m*VQ; H=(2*n+2)*VP+h*VQ; pi=k*VQ; O2=C+H/4; M=(12*C+H)/10
    print(C,H,pi,O2,M); assert (C,H,pi,O2,M)==(35,60,15,50,48)
    ratio=VP/VQ if m==4 else None; print('alkane ratio',ratio)
assert out[0][4]/out[0][5]==F(1,3) and out[1][4]/out[1][5]==1
for e in ['C4H4 + 3H2 -> C4H10','C4H6 + 2H2 -> C4H10','C3H8 + 5O2 -> 3CO2 + 4H2O','C4H4 + 5O2 -> 4CO2 + 2H2O']:
    assert check_equation(e) is None

print('全検証OK:', ITEM)
