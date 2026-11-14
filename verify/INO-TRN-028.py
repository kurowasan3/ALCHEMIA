# ---- common: check every complete \ce{} equation in an item ----
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
def ce_blocks(s):
    out=[]; i=0
    while True:
        j=s.find(r'\ce{', i)
        if j<0: return out
        k=j+4; depth=1
        while depth:
            if s[k]=='{': depth+=1
            elif s[k]=='}': depth-=1
            k+=1
        out.append(s[j+4:k-1]); i=k
def check_item_equations(path):
    it=json.load(open(path, encoding='utf-8'))
    n=0
    for key in ('problem','answer','solution'):
        for b in ce_blocks(it[key]):
            arrows=len(re.findall(r'->|<=>', b))
            if arrows!=1 or '\\boxed' in b or '\\text' in b or '\\mathrm' in b: continue
            b2=re.sub(r'\((s|l|g|aq)\)','',b)
            r=check_equation(b2.replace('<=>','->'))
            print(('OK ' if r is None else 'NG ')+b, '' if r is None else r)
            assert r is None, b
            n+=1
    print('equations checked:', n)
    return it
from fractions import Fraction as F
from itertools import combinations
it=check_item_equations('/home/claude/chem3/outh/items/INO-TRN-028/item.json')
c=F(2,100)
n2=5*c*F(1020,100000); nT=5*c*F(5102,100000); n3=nT-n2
nO=(2*n2+3*n3)/2
mass=56*nT+16*nO
print(float(n2),float(n3),float(nO),float(mass)); assert abs(mass-F(4,10))<F(1,10000)
# oxides: (Fe2+, Fe3+, O, M)
ox={'FeO':(1,0,1,72),'Fe3O4':(1,2,4,232),'Fe2O3':(0,2,3,160)}
sols=[]
for k in ox:  # single
    a,b,o,M=ox[k]
    if a*n3==b*n2 and (a+b)>0: sols.append((k,))
for k1,k2 in combinations(ox,2):
    a1,b1,o1,M1=ox[k1]; a2,b2,o2,M2=ox[k2]
    det=a1*b2-a2*b1
    if det==0: continue
    x=(n2*b2-a2*n3)/det; y=(a1*n3-b1*n2)/det
    if x>0 and y>0:
        sols.append((k1,float(x*M1),k2,float(y*M2)))
        assert abs(x*o1+y*o2-nO)<F(1,10**9)
print(sols)
names=sorted(tuple(sorted([s[0],s[2]])) for s in sols)
assert names==[('FeO','Fe2O3'),('Fe2O3','Fe3O4')] or names==sorted([('Fe2O3','FeO'),('Fe2O3','Fe3O4')]), names
d={ (s[0],s[2]):(round(s[1],4),round(s[3],4)) for s in sols}
print(d)
heat=nT/2*160; print('air heated', float(heat)); assert round(float(heat),3)==0.408
print('ALL OK')
