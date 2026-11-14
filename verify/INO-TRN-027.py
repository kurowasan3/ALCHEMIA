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
it=check_item_equations('/home/claude/chem3/out/items/INO-TRN-027/item.json')
cuo=F(240,100)/80; h2o=F(18,100)/18; co2=F(88,100)/44
print(cuo,h2o,co2); assert co2/h2o==2 and co2+h2o==cuo
M=2*(64+12+48)+(64+34); assert M==346 and M*h2o==F(346,100)
print('Cu g', float(cuo*64)); assert cuo*64==F(192,100)
print('ALL OK')
