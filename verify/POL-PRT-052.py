ITEM = '/home/claude/chemx/kyodai/out/items/POL-PRT-052/item.json'
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

from rdkit import Chem
chk('rise', 0.54/3.6, 0.15, 2)
n = 1.32e4/110; assert round(n) == 120
h = 0.60*n; assert round(h) == 72
chk('長さ', h*0.15, 10.8, 3); assert round(h/3.6) == 20
assert round(h) - 4 == 68
chk('比', 0.35/0.15, 2.3, 2)
dip = Chem.AddHs(Chem.MolFromSmiles('CC(N)C(=O)N1CCCC1C(=O)O'))  # Ala-Pro
amideN = [a for a in dip.GetAtoms() if a.GetSymbol()=='N' and a.IsInRing()][0]
assert amideN.GetTotalNumHs() == 0 and not any(nb.GetSymbol()=='H' for nb in amideN.GetNeighbors())

print('ALL OK')
