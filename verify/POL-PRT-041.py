ITEM = '/home/claude/chem5/out/items/POL-PRT-041/item.json'
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
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
A = {'H':1,'C':12,'N':14,'O':16}
def fm(c,h,n,o): return c*12+h*1+n*14+o*16
res = {'Gly':(2,3,1,1),'Ala':(3,5,1,1),'Lys':(6,12,2,1),'Phe':(9,9,1,1),'Arg':(6,12,4,1)}
aa = {'Gly':'NCC(=O)O','Ala':'NC(C)C(=O)O','Lys':'NCCCCC(N)C(=O)O','Phe':'NC(Cc1ccccc1)C(=O)O','Arg':'NC(=N)NCCCC(N)C(=O)O'}
import re as _re
for k,s in aa.items():
    f = CalcMolFormula(Chem.MolFromSmiles(s))
    d = {e:int(n or 1) for e,n in _re.findall(r'([A-Z][a-z]?)(\d*)', f)}
    c,h,n,o = res[k]
    assert (d['C'],d['H']-2,d['N'],d['O']-1) == (c,h,n,o), (k,f)
print('残基の組成式 OK')
M = {k:fm(*v) for k,v in res.items()}; print(M)
assert sig(100*14/M['Gly'],3)==24.6 and sig(100*14/M['Phe'],3)==9.52
frac = {'Gly':R(40),'Ala':R(30),'Lys':R(10),'Phe':R(20)}
mass = sum(frac[k]*M[k] for k in frac); nN = sum(frac[k]*res[k][2] for k in frac)
Np = float(nN*14/mass*100); fac = 100/Np
print('N%', Np, '係数', fac)
assert sig(Np,3)==17.8 and sig(fac,3)==5.60
N = 0.170*14/17
p1 = N*6.25/10.0*100; p2 = N*sig(fac,3)/10.0*100; p2b = N*fac/10*100
print(p1,p2,p2b); assert sig(p1,2)==8.8 and sig(p2,2)==7.8 and sig(p2b,2)==7.8
print('Arg 残基 N%', 100*4*14/M['Arg']); assert sig(100*4*14/M['Arg'],3)==35.9
print('ALL OK')
