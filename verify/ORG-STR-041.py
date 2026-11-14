import sys, re, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def ce_blocks(s):
    out = []; i = 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq): return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq): return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1: return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_all_equations(item_path):
    d = json.load(open(item_path, encoding='utf-8'))
    n = 0
    for k in ('problem', 'answer', 'solution'):
        for b in ce_blocks(d[k]):
            if should_check(b):
                r = check_equation(clean_ce(b)); n += 1
                assert r is None, (b, r)
    print('equations checked:', n)
    return d
from sympy import Rational as Q, nsimplify
d = check_all_equations('/home/claude/chem5/out/items/ORG-STR-041/item.json')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
# 組成式
C,Hh=79.2,5.7; O=100-C-Hh
r=[C/12,Hh/1.0,O/16]; m=min(r); print('ratio',[x/m for x in r])
assert all(abs(a-b)<0.1 for a,b in zip([x/m for x in r],[7,6,1]))
# 想定物質（安息香酸ベンジル）と組成の一致
mw=14*12+12*1.0+2*16; print('MW C14H12O2',mw)
print('calc %C',14*12/mw*100,'%H',12/mw*100)
assert CalcMolFormula(Chem.MolFromSmiles('O=C(OCc1ccccc1)c1ccccc1'))=='C14H12O2'
M=40*0.0636/(1.20*0.0100); print('M',M); assert abs(M-212)<0.5
dtb=5.12*(0.0636/212)/0.0100; print('benzene dT',dtb); assert round(dtb,2)==0.15
# 誤差範囲
cam=(40*0.0636/0.0100/1.30, 40*0.0636/0.0100/1.10); print('camphor M range',cam)
ben=(5.12*0.0636/0.0100/(dtb+0.10), 5.12*0.0636/0.0100/(dtb-0.10)); print('benzene M range',ben)
nc=[n for n in range(1,10) if cam[0]<=106*n<=cam[1]]; nb=[n for n in range(1,10) if ben[0]<=106*n<=ben[1]]
print(nc,nb); assert nc==[2] and nb==[2,3,4,5]
print('rel err camphor',0.10/1.20,'benzene',0.10/dtb)
print('OK')
