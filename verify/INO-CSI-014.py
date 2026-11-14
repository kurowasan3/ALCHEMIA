# 検証スクリプト INO-CSI-014  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/INO-CSI-014/item.json'
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
Na,Al,Si,O,H,Ca = 23.0,27.0,28.1,16.0,1.0,40.1
X = 8*Na+8*Al+40*Si+96*O+24*18.0; assert abs(X-3492) < 1e-6
assert 2*(8+40) == 96 and 2*(12+12) == 48
assert sig(432/X*100,3) == 12.4
assert sig(4/X*Ca*1e3,2) == 46
Y = 12*Na+12*Al+12*Si+48*O+27*18.0; assert abs(Y-2191.2) < 1e-6
capY = 6/Y*Ca*1e3; assert sig(capY,2) == 110
assert sig(80/capY,2) == 0.73
print('ALL OK')
