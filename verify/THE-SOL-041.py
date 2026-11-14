# 検証スクリプト THE-SOL-041  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-SOL-041/item.json'
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

RT = 8.31e3*300
pi_salt = 0.300*RT; print('NaCl', pi_salt); assert sig(pi_salt,2) == 7.5e5
col = 76.0*pi_salt/1.013e5*13.6/1.03/100; print('液柱 m', col); assert round(col) == 74
p = 1.013e5*(7.48*1.03/13.6)/76.0; print('pi', p); assert sig(p,2) == 7.6e2 and round(p) == 755
c = p/RT; M = 20.0/c; print('M', M); assert sig(M,2) == 6.6e4
ratio = p/(pi_salt+p); print('割合', ratio); assert 0.0009 < ratio < 0.0011

print('全検証OK:', ITEM)
