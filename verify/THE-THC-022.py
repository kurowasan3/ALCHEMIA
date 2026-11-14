# 検証スクリプト THE-THC-022  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-THC-022/item.json'
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

from sympy import symbols, solve, Rational as Q
HH, OO, C, CH4 = symbols('HH OO C CH4')
s = solve([Q(1,2)*946 + Q(3,2)*HH - 1170 + 46,
           HH + OO/2 - 928 + 242,
           C + OO - 1608 + 394,
           C + 2*HH - CH4 + 75], [HH, OO, C, CH4]); print(s)
assert s[HH] == 434 and s[OO] == 504 and s[C] == 710 and s[CH4] == 1653
assert round(s[CH4]/4) == 413
co = s[C] + s[OO]/2 - 1074; print('CO', co); assert co == -112
for e in ['1/2N2 + 3/2H2 -> NH3','H2 + 1/2O2 -> H2O','C + O2 -> CO2','C + 2H2 -> CH4','C + 1/2O2 -> CO']: assert check_equation(e) is None

print('全検証OK:', ITEM)
