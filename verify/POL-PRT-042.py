ITEM = '/home/claude/chem5/out/items/POL-PRT-042/item.json'
# 検証スクリプト（グループ g09）  python3 -I この.py
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert [int(m.group(2)) for m in re.finditer(r'(^|\s)\((\d+)\)', it['answer'])] == list(range(1, np_+1))

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from sympy import symbols, solve, Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
bases = {'A':'Nc1ncnc2[nH]cnc12','G':'Nc1nc2[nH]cnc2c(=O)[nH]1','C':'Nc1cc[nH]c(=O)n1','T':'Cc1c[nH]c(=O)[nH]c1=O','U':'O=c1cc[nH]c(=O)[nH]1'}
want = {'A':'C5H5N5','G':'C5H5N5O','C':'C4H5N3O','T':'C5H6N2O2','U':'C4H4N2O2'}
nN = {}
for k,s in bases.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m)==want[k]
    nN[k] = sum(1 for a in m.GetAtoms() if a.GetSymbol()=='N')
print(nN)
assert nN['A']+nN['T']==7 and nN['G']+nN['C']==8
f = symbols('f')
expr = (nN['A']+nN['T'])*(1-f)/2 + (nN['G']+nN['C'])*f/2
assert (expr - (7+f)/2).simplify()==0
fv = solve(expr-R('3.70'),f)[0]; print('f=',fv); assert fv==R(2,5)
for comp in [dict(A=.25,G=.25,C=0,U=.5), dict(A=.2,G=.2,C=.3,U=.3)]:
    v = sum(comp[k]*nN[k] for k in comp); print(comp, v); assert abs(v-3.5)<1e-12 and abs(sum(comp.values())-1)<1e-12
print('ALL OK')
