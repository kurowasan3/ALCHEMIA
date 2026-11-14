ITEM = '/home/claude/chemx/kyodai/out/items/POL-SUG-048/item.json'
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
pat = Chem.MolFromSmarts('[OX2H]')
cands = {
 'ア': ('OC1OC(CO)C(OC2OC(CO)C(O)C(OC3OC(CO)C(O)C(O)C3O)C2O)C(O)C1O', 'C18H32O16', 11, (3,0,0)),
 'イ': ('OC1OC(CO)C(OC2OC(CO)C(O)C(O)C2OC3OC(C)C(O)C(O)C3O)C(O)C1O', 'C18H32O15', 10, (2,1,0)),
 'ウ': ('OC1OC(CO)C(OC2OC(CO)C(O)C(OC3OC(CO)C(O)C(O)C3NC(C)=O)C2O)C(O)C1O', 'C20H35NO16', 10, (2,0,1)),
 'エ': ('OC1OC(CO)C(OC2OC(CO)C(O)C(OC3OC(C)C(O)C(O)C3O)C2NC(C)=O)C(O)C1O', 'C20H35NO15', 9, (1,1,1)),
}
ok = []
for k, (smi, f, noh, (a,b,c)) in cands.items():
    m = Chem.MolFromSmiles(smi)
    assert CalcMolFormula(m) == f, (k, CalcMolFormula(m))
    n = len(m.GetSubstructMatches(pat)); assert n == noh == 3*a+2*b+2*c+2, (k, n)
    M = round(sum({'C':12,'H':1,'O':16,'N':14}[at.GetSymbol()] for at in Chem.AddHs(m).GetAtoms()))
    need = 0.488/M*noh
    print(k, f, M, noh, need)
    if abs(need-1.42/142) < 2e-5: ok.append(k)
assert ok == ['イ'], ok
malt = Chem.MolFromSmiles('OC1OC(CO)C(OC2OC(CO)C(O)C(O)C2O)C(O)C1O'); assert len(malt.GetSubstructMatches(pat)) == 8

print('ALL OK')
