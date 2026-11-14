# 検証スクリプト THE-THC-023  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-THC-023/item.json'
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
    return re.sub(r'\{(\d*[+-])\}', r'\1', eq)
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'<=>|->|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                e = check_equation(clean_ce(b)); n += 1
                assert e is None, (b, e)
    print('反応式', n, '本 つり合いOK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

Cp = 100*4.18
nMg = 0.243/24.3; nMgO = 0.403/40.3; assert abs(nMg-0.01)<1e-12 and abs(nMgO-0.01)<1e-12
dA = -Cp*11.2/nMg/1000; dB = -Cp*3.62/nMgO/1000; print(dA, dB)
assert sig(dA,3) == -468 and sig(dB,3) == -151
dC = dA - dB - 286; print(dC); assert sig(dC,3) == -603
nHCl = 0.0500; n3 = min(0.729/24.3, nHCl/2); assert abs(n3-0.025)<1e-12
dT = n3*(-dA)*1000/Cp; print(dT); assert sig(dT,3) == 28.0
left = (0.729/24.3 - n3)*24.3; print(left); assert abs(left-0.1215)<1e-9 and sig(left,2) == 0.12
for e in ['Mg + 2HCl -> MgCl2 + H2','MgO + 2HCl -> MgCl2 + H2O','H2 + 1/2O2 -> H2O','Mg + 1/2O2 -> MgO']: assert check_equation(e) is None

print('全検証OK:', ITEM)
