ITEM = '/home/claude/chemx/kyodai/out/items/ORG-STR-059/item.json'
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
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Descriptors
for smi, f, m in [('CC(O)c1ccccc1','C8H10O',122),('CC(=O)OC=C','C4H6O2',86),('CC(=O)OC(C)c1ccccc1','C10H12O2',164)]:
    mol = Chem.MolFromSmiles(smi); assert CalcMolFormula(mol) == f; assert abs(Descriptors.MolWt(mol)-m) < 0.5
assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles('CC(O)c1ccccc1'), includeUnassigned=True)) == 1
n = 12.2/122; v = 3.44/86
chk('エステル', v*164, 6.56, 3); chk('R残', (n/2-v)*122, 1.22, 3); chk('S残', n/2*122, 6.10, 3)
chk('ee', (n/2-(n/2-v))/(n/2+(n/2-v))*100, 66.7, 3)
chk('最小', n/2*86, 4.30, 3); chk('エステル2', n/2*164, 8.20, 3)
assert 6.00/86 > n/2

print('ALL OK')
