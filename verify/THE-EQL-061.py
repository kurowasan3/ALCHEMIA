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
from sympy import symbols, solve, Rational as R
it=check_item_equations('/home/claude/chem3/outh/items/THE-EQL-061/item.json')
K1,K2=R(3,2),R(1,2)
def phase(co,co2):
    r=R(co2,co)
    return 'Fe3O4' if r>K1 else ('FeO' if r>K2 else 'Fe')
print([phase(60,40),phase(75,25),phase(30,70)]); assert [phase(60,40),phase(75,25),phase(30,70)]==['FeO','Fe','Fe3O4']
x=symbols('x'); s=solve(x/(2-x)-K2,x)[0]; print('Fe',s,float(s)); assert s==R(2,3)
assert 1+1/K2==3
assert 4+4/K2==12
n=symbols('n'); nlow=solve(3/(n-3)-K2,n)[0]; nup=solve(4/(n-4)-K1,n)[0]
print(nlow, float(nup)); assert nlow==9 and nup<nlow
r_out=R(4,5); assert r_out<=K1
print('ALL OK')
