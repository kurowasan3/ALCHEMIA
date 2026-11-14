# 検証スクリプト THE-SOL-077  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-SOL-077/item.json'
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', b): continue  # 一般式（validate と同じ除外）
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式', n, '本 OK')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure(); check_all_equations()
from math import sqrt, log10
g = lambda z, C: 10**(-0.5*z**2*sqrt(z**2*C))
gk, gm = g(1, 1e-2), g(2, 1e-3)
assert (sig(gk,3), sig(gm,3)) == (0.891, 0.747)
assert (sig(2*gk,3), sig(1+gk,3), sig(2*gm,3), sig(1+gm,3)) == (1.78, 1.89, 1.49, 1.75)
vals = {z: 2*g(z, 4.0e-3) for z in (1,2,3)}; print(vals)
assert [z for z in vals if abs(vals[z]-1.12) < 0.02] == [2]
assert sig(1+g(2,4.0e-3), 2) == 1.6
th = {z: (log10(0.90)/(-0.5*z**3))**2 for z in (1,2)}; print(th)
assert sig(th[1],2) == 8.4e-3 and sig(th[2],2) == 1.3e-4 and abs(th[1]/th[2]-64) < 1e-9
gI = lambda z, I: 10**(-0.5*z**2*sqrt(I))
assert sig(gI(1,1e-2),2) == 0.89 and sig(gI(2,1e-2),2) == 0.63
print('ALL OK')
