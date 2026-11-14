ITEM = '/home/claude/chemx/kyodai/out/items/THE-BND-074/item.json'
# 検証スクリプト（京大 g01）  python3 -I この.py
import sys, json, re, math
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)) or v == expect, (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

# 小問の数と答えの数の一致
nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

d = 4*144/(6.0e23*0.182e-21); chk('密度', d, 5.27, 3)
r = 1.73/4*0.567; chk('距離', r, 0.245, 3)
chk('密度比(厳密)', (2/3**0.5)**3, 1.5, 2)
chk('密度比(√3=1.73)', (2/1.73)**3, 1.5, 2)
chk('密度比(数値)', (0.567/(2*r))**3, 1.5, 2)

print('ALL OK')
