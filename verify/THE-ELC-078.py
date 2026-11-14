ITEM = '/home/claude/chemx/x1/out/items/THE-ELC-078/item.json'
# 検証スクリプト（x1 g04）  python3 -I この.py
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
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

F=9.65e4; e=1.00*1930/F; chk('e',e,2.00e-2,3)
chk('Cu',e/2*63.5,0.635,3)
O2=e/4; dis=O2-4.78e-3; chk('dis',dis,2.20e-4,3)
p=4.78e-3/(4.78e-3+1.00e-2); chk('p',p,0.323,3)
s=dis/(0.500*p); chk('s',s,1.36e-3,3)
d4=s*1.00*0.500; chk('d4',d4,6.80e-4,3)
chk('比',d4/dis,3.1,2)

print('ALL OK')
