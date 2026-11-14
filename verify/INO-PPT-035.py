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
it=check_item_equations('/home/claude/chem3/out/items/INO-PPT-035/item.json')
assert 137+32+64==233 and 137+12+48==197
so4=F(233,100)/233; co3=(F(430,100)-F(233,100))/197
print(so4,co3); assert so4==co3==F(1,100)
print('CO2 L', float(co3*F(224,10))); assert co3*F(224,10)==F(224,1000)
import json
ions=json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json',encoding='utf-8'))
s=json.dumps(ions,ensure_ascii=False)
for f,c in [('BaSO4','白'),('BaCO3','白'),('BaCrO4','黄')]:
    m=re.search(r'"f": "%s", "color": "([^"]+)"'%f, s); print(f, m.group(1) if m else None)
    assert m and m.group(1)==c
print('ALL OK')
