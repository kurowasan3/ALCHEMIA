ITEM = '/home/claude/chemx/kyodai/out/items/POL-FIB-019/item.json'
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
for e in ['HOOC(CH2)4COOH + H2N(CH2)6NH2 -> HOOC(CH2)4CONH(CH2)6NH2 + H2O',
          'ClOC(CH2)4COCl + H2N(CH2)6NH2 -> ClOC(CH2)4CONH(CH2)6NH2 + HCl']:
    assert check_equation(e) is None, e
AW = {'C':12,'H':1.0,'N':14,'O':16,'Cl':35.5}
def mw(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); return sum(AW[a.GetSymbol()] for a in m.GetAtoms())
assert mw('OC(=O)CCCCC(=O)O') == 146 and mw('ClC(=O)CCCCC(=O)Cl') == 183 and mw('COC(=O)CCCCC(=O)OC') == 174 and mw('NCCCCCCN') == 116
def unit_formula(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); c = {}
    for a in m.GetAtoms():
        if a.GetSymbol() != '*': c[a.GetSymbol()] = c.get(a.GetSymbol(), 0) + 1
    return c
u66 = unit_formula('*C(=O)CCCCC(=O)NCCCCCCN*'); assert u66 == {'C':12,'H':22,'N':2,'O':2}, u66
ut = unit_formula('*C(=O)C(O)C(O)C(=O)NCCCCCCN*'); assert ut == {'C':10,'H':18,'N':2,'O':4}, ut
assert sum(AW[k]*v for k,v in u66.items()) == 226 and sum(AW[k]*v for k,v in ut.items()) == 230
n = 1000/226
chk('acid', n*146, 646, 3); chk('H2O', n*36, 159, 3)
chk('chloride', n*183, 810, 3); chk('HCl', n*73, 323, 3)
chk('ester', n*174, 770, 3); chk('MeOH', n*64, 283, 3)
chk('NaOH', 2*n, 8.85, 3)
assert abs(n*146 + n*116 - (1000 + n*36)) < 1e-6

print('ALL OK')
