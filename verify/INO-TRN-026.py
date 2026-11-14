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
it=check_item_equations('/home/claude/chem3/out/items/INO-TRN-026/item.json')
M=56+2*52+4*16; assert M==224
# electrons lost per FeCr2O4
e=(3-2)+2*(6-3); assert e==7
n=F(1000*448,1000)/M; assert n==2
o2=n*F(7,4); print('O2 L',float(o2*F(224,10))); assert o2*F(224,10)==F(784,10)
cr=2*n*52; al=2*n*27; print(cr,al); assert cr==208 and al==108
print('ALL OK')
