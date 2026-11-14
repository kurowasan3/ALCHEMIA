ITEM = '/home/claude/chem6/out/items/INO-NIT-018/item.json'
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
# 工程X
hno3 = F(2,3); o2 = F(5,4)+F(1,2)
assert hno3 == F(2,3) and o2/hno3 == F(21,8)
print('X: HNO3', hno3, 'O2/HNO3', o2/hno3)
# 全体式 = (①+②×3+③×2)/4  係数ベクトル {NH3,O2,NO,H2O,NO2,HNO3}
import collections
def add(*terms):
    c=collections.Counter()
    for k,eq in terms:
        for s,v in eq.items(): c[s]+=k*v
    return {s:v for s,v in c.items() if v}
e1={'NH3':-4,'O2':-5,'NO':4,'H2O':6}; e2={'NO':-2,'O2':-1,'NO2':2}; e3={'NO2':-3,'H2O':-1,'HNO3':2,'NO':1}
tot=add((1,e1),(3,e2),(2,e3)); print(tot)
assert tot=={'NH3':-4,'O2':-8,'H2O':4,'HNO3':4}
for e in ['4NH3 + 5O2 -> 4NO + 6H2O','2NO + O2 -> 2NO2','3NO2 + H2O -> 2HNO3 + NO','NH3 + 2O2 -> HNO3 + H2O','4NO + 3O2 + 2H2O -> 4HNO3']:
    assert check_equation(e) is None
n = 3.4e6/17.0; mY = n*63.0/0.63/1e6; mX = mY*2/3
print('Y', mY, 'X', mX); assert sig(mY,2)==20 and sig(mX,2)==13

print('全検証OK:', ITEM)
