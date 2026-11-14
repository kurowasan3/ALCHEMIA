# 検証スクリプト THE-THC-021  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-THC-021/item.json'
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

dc = (-394 + 2*(-286)) - (-74); assert dc == -892
dh = 54 + dc; assert dh == -838
M = 12 + 4*1.0 + 6*(2*1.0+16); assert M == 124
q1 = 1000/M*838; q2 = 1000/16.0*892; print(q1, q2); assert sig(q1,3) == 6.76e3 and sig(q2,3) == 5.58e4
n = 910/M; v = n*22.4; print(n, v); assert sig(v,3) == 164
x = 54/892; print(x); assert sig(100*x,2) == 6.1
assert check_equation('CH4*6H2O + 2O2 -> CO2 + 8H2O') is None

print('全検証OK:', ITEM)
