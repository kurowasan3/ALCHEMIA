ITEM = '/home/claude/chemx/x1/out/items/POL-RES-029/item.json'
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

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
AW={'H':1.0,'C':12,'O':16,'Cl':35.5}
def mw(m):
    m=Chem.AddHs(m); return sum(AW[a.GetSymbol()] for a in m.GetAtoms())
G='C1OC1C'; Ar='c1ccc(cc1)C(C)(C)c1ccc(cc1)'
for n in range(4):
    smi='C1OC1COc1ccc(cc1)C(C)(C)c1ccc(cc1)O' + 'CC(O)COc1ccc(cc1)C(C)(C)c1ccc(cc1)O'*n + 'CC1CO1'
    m=Chem.MolFromSmiles(smi)
    assert abs(mw(m)-(284*n+340))<1e-9, (n,mw(m))
    if n==0: assert CalcMolFormula(m)=='C21H24O4'
print('分子量 284n+340 OK')
assert abs(mw(Chem.MolFromSmiles('CC(C)(c1ccc(O)cc1)c1ccc(O)cc1'))-228)<1e-9
assert abs(mw(Chem.MolFromSmiles('ClCC1CO1'))-92.5)<1e-9
chk('n',(2*525-340)/284,2.5,2)
from fractions import Fraction
r=Fraction(35,10)/Fraction(45,10); assert (r.numerator,r.denominator)==(7,9)

print('ALL OK')
