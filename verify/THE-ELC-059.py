# 検証スクリプト（グループ g08）  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-ELC-059/item.json'
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
print('反応式のつり合い OK：', n_eq, '式')
# 小問数と答えの数の一致
m_en = re.search(r'\\begin\{enumerate\}(.*)\\end\{enumerate\}', it['problem'], re.S)
nq = m_en.group(1).count('\\item') if m_en else 0
assert '\\begin{itemize}' not in it['problem']
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問数', nq, '答えの数', na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


from sympy import Rational as R
e = 4*R(25,1000)
eH2, eCO, eCH4 = 2*R(15,1000), 2*R(10,1000), 8*R(25,10000)
fe = [x/e*100 for x in (eH2,eCO,eCH4)]; print(e, fe); assert fe == [30,20,20]
eX = e - eH2 - eCO - eCH4; cX = R(275,10000) - R(10,1000) - R(25,10000)
print(eX, cX, eX/cX); assert eX == R(3,100) and cX == R(15,1000) and eX/cX == 2
cand = {'HCOOH':2, 'HCHO':4, 'CH3OH':6}
assert [k for k,v in cand.items() if v == eX/cX] == ['HCOOH']
for eq in ['CO2 + 2H+ + 2e- -> HCOOH', 'CO2 + 4H+ + 4e- -> HCHO + H2O', 'CO2 + 6H+ + 6e- -> CH3OH + H2O', '2H+ + 2e- -> H2']:
    assert check_equation(eq) is None, eq
tot = 2/(eX/e); print(float(tot), float(tot*96500), float(tot/4))
assert sig(float(tot*96500),2) == 640000 and sig(float(tot/4),2) == 1.7

print('ALL OK')
