ITEM = '/home/claude/chemx/tokodai/out/items/POL-RES-026/item.json'
# 検証スクリプト（東工大 g04）  python3 -I この.py
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
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

from fractions import Fraction as F
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
mass = {'C':12,'H':1,'N':14,'O':16}
def fw(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    return sum(mass[a.GetSymbol()] for a in m.GetAtoms())
# 単位を末端 H で近似せず，主鎖の 2 炭素を含む単位として：CH2-CH(Ar) = C8H8 相当（エチルベンゼン − 2H）
u0 = fw('CCc1ccccc1')-2; u1 = fw('CCc1ccc([N+](=O)[O-])cc1')-2; u2 = fw('CCc1ccc(N)cc1')-2; u3 = fw('CCc1ccc(NC(C)=O)cc1')-2
print(u0,u1,u2,u3)
assert (u0,u1,u2,u3) == (104,149,119,161)
n = F(104,10)/104
x = (F(131,10)-F(104,10))/(u1-u0)
chk('(2)', float(x/n*100), 60, 2)
chk('(3)', float((n-x)*u0 + x*u3), 13.8, 3)
chk('Y', float((n-x)*u0 + x*u2), 11.3, 3)

print('ALL OK')
