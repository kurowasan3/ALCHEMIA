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
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from sympy import Rational as Q, symbols, expand
d = check_all_equations('/home/claude/chem5/outh/items/ORG-STR-047/item.json')
def stereo(smi):
    m=Chem.MolFromSmiles(smi)
    iso={Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m,StereoEnumerationOptions(unique=True))}
    def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@','!').replace('@','@@').replace('!','@')))
    achiral=[s for s in iso if mirror(s)==s]
    return iso,achiral
Z='CCC(C)OC(=O)c1ccccc1C(=O)OC(C)CC'
iso,ach=stereo(Z); print('Z',len(iso),'achiral',len(ach)); assert len(iso)==3 and len(ach)==1
print(CalcMolFormula(Chem.MolFromSmiles(Z)))
Z4='CCC(C)OC(=O)c1ccc(C)cc1C(=O)OC(C)CC'
iso4,ach4=stereo(Z4); print('Z4',len(iso4),'achiral',len(ach4)); assert len(iso4)==4 and len(ach4)==0
# 確率
p=symbols('p'); q=1-p
for pv in (Q(1,2),Q(3,4)):
    RR,SS,RS=(p**2).subs(p,pv),(q**2).subs(p,pv),(2*p*q).subs(p,pv)
    print(pv,RR,SS,RS, 'sum',RR+SS+RS)
pv=Q(3,4); RR,SS,RS=pv**2,(1-pv)**2,2*pv*(1-pv)
assert (RR,SS,RS)==(Q(9,16),Q(1,16),Q(6,16)) and RS*100==Q(75,2)
ratio=RR/SS; print('R:S after',ratio); assert ratio==9
ee0=(Q(3,4)-Q(1,4)); ee1=(RR-SS)/(RR+SS); print('ee',ee0,ee1); assert ee0==Q(1,2) and ee1==Q(4,5)
fracS_meso=(RS*1)/(RS*1+2*SS)  # meso は S を1個, (S,S) は2個; print('S in meso',fracS_meso); assert fracS_meso==Q(3,4)
print('OK')
