ITEM = '/home/claude/chemx/x1/out/items/INO-PPT-049/item.json'
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

for e in ['Cu^2+ + 2NH3 + 2H2O -> Cu(OH)2 + 2NH4+','Cu(OH)2 + 4NH3 -> [Cu(NH3)4]^2+ + 2OH-']:
    assert check_equation(e) is None, e
# 平面四角形 4 位置（0..3 を一周），配位子は隣り合う 2 位置を占める
import itertools
def configs(ends):
    res=set()
    for o in itertools.product(range(2), repeat=2):
        # 配位子1: 位置0,1  配位子2: 位置2,3
        a=[None]*4
        a[0],a[1] = (ends[0],ends[1]) if o[0]==0 else (ends[1],ends[0])
        a[2],a[3] = (ends[0],ends[1]) if o[1]==0 else (ends[1],ends[0])
        res.add(tuple(a))
    def canon(a):
        ims=[]
        for r in range(4):
            b=a[r:]+a[:r]; ims += [b, b[::-1]]
        return min(ims)
    return {canon(a) for a in res}
assert len(configs(('N','O')))==2 and len(configs(('N','N')))==1
print('幾何異性体の数OK')

print('ALL OK')
