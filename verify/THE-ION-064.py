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
from sympy import symbols, solve, sqrt, log, N
it=check_item_equations('/home/claude/chem3/outh/items/THE-ION-064/item.json')
Ka,Ksp,Kw=1e-5,1e-33,1e-14
import math
def exact(C): return (-Ka+math.sqrt(Ka*Ka+4*Ka*C))/2
x1=exact(0.10); print('x1',x1,-math.log10(x1)); assert round(-math.log10(math.sqrt(Ka*0.1)),1)==3.0 and abs(x1-1e-3)/1e-3<0.01
Q1=(0.1-x1)*(Kw/x1)**3; print('Q1',Q1); assert Q1<Ksp
pH_on=14+math.log10((Ksp/0.1)**(1/3)); print('pH onset',pH_on); assert round(pH_on,1)==3.3
x3=exact(1e-4); approx=math.sqrt(Ka*1e-4); print('x3',x3,'approx',approx,'pH',-math.log10(x3), 'alpha',x3/1e-4)
assert round(-math.log10(x3),1)==4.6 and abs((-1+6.40)/2*1e-5-x3)/x3<0.01
Q3=(1e-4-x3)*(Kw/x3)**3; print('Q3',Q3); assert Q3>Ksp and round(Q3/1e-33,1)==3.7
print('5-0.43 =',5-0.43)
print('ALL OK')
