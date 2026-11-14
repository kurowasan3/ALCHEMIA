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
it=check_item_equations('/home/claude/chem3/out/items/INO-TRN-024/item.json')
Fe,O=56,16
nFe=F(168,100)/Fe; nO=(F(232,100)-F(168,100))/O
print(nFe,nO); assert nFe/nO==F(3,4)
fr={'FeO':F(56,72),'Fe3O4':F(168,232),'Fe2O3':F(112,160)}
obs=F(168,232)
assert [k for k,v in fr.items() if v==obs]==['Fe3O4']
m=nFe/2*160; print('Fe2O3 mass',float(m)); assert m==F(240,100)
ratio=F(232,72); print('ratio',float(ratio)); assert round(float(ratio),1)==3.2
print('ALL OK')
