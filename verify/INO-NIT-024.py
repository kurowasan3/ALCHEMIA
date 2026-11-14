# 検証スクリプト INO-NIT-024  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/INO-NIT-024/item.json'
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
sols = [(a,b,c,d) for a in range(1,9) for b in range(1,9) for c in range(1,9) for d in range(1,13)
        if check_equation(f'{a}NO + {b}NH3 + O2 -> {c}N2 + {d}H2O') is None]
print(sols); assert sols == [(4,4,4,6)]
assert check_equation('6NO + 4NH3 -> 5N2 + 6H2O') is None
nh3 = 2.0*4/4 + 0.50*8/6; urea = nh3/2*60.0; print(urea); assert sig(urea,2) == 80
assert abs((8+4)/8 - 1.5) < 1e-12
print('ALL OK')
