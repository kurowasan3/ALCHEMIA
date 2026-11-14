# 検証スクリプト（グループ g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_all_equations(it, extra=()):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|…', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]', b): continue
            if not (' + ' in b or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', b)): continue
            e = clean_ce(b)
            r = check_equation(e)
            assert r is None, (b, r)
            print('式OK:', b); n += 1
    for e in extra:
        r = check_equation(e); assert r is None, (e, r); print('式OK(追加):', e); n += 1
    return n

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def check_parts(it):
    k = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, k + 1)], (k, nums)
    print('小問数OK:', k)
ITEM = '/home/claude/chemk/out/items/THE-MOL-048/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); assert check_all_equations(it) >= 1
def V(V0, a, v):
    no_react = min(v, 2*a); o2_react = no_react/2
    return V0 + v - no_react - o2_react
tab = [(0,500),(40,480),(80,460),(120,440),(160,420),(200,460),(240,500)]
ok = [a/2 for a in range(0, 1001) if all(abs(V(500, a/2, v) - W) < 1e-9 for v, W in tab)]
print('試料1 a =', ok); assert ok == [80.0]; print('割合', 80/500*100)
assert abs(V(500, 105, 300) - 485) < 1e-9 and 500 - 300/2 == 350
a2 = (500 + 300 - 485)/3; assert a2 == 105 and 2*a2 <= 300; print('試料2', a2/500*100)
good = [a/2 for a in range(0, 1001) if abs(V(500, a/2, 150) - 425) < 1e-9]
print('試料3 a の範囲', min(good), max(good)); assert min(good) == 75.0 and max(good) == 500.0
print('全検証OK:', ITEM)
