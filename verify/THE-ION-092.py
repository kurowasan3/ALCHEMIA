# 検証スクリプト THE-ION-092  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/THE-ION-092/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


from math import log10
def start(Ksp, n, c=1e-2):
    return (log10(Ksp) + 14*n - log10(c))/n
r = [round(start(1e-11,2),2), round(start(1e-15,2),2), round(start(1e-38,3),2)]
print(r); assert r == [9.5, 7.5, 2.0]
n = (-3 - (-6))/(5.0-4.0); assert n == 3
logK = -3 - 14*3 + 3*4.0; print(logK); assert abs(logK + 33) < 1e-12
print('OK')
