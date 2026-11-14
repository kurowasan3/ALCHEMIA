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
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
d = check_all_equations('/home/claude/chem5/out/items/ORG-CBX-032/item.json')
def comp(C,H):
    O=100-C-H; r=[C/12,H/1.0,O/16]; m=min(r); return [x/m for x in r]
print('P ratio',comp(40.7,5.1),'Q ratio',comp(41.4,3.4))
eqP=0.177/(0.100*30.0/1000); eqQ=0.145/(0.100*25.0/1000); print('eq',eqP,eqQ)
assert abs(eqP-59.0)<0.05 and abs(eqQ-58.0)<0.05
# P: (C2H3O2)_k, k COOH; H偶数・不飽和度>=k
okP=[k for k in range(1,9) if (3*k)%2==0 and (2*2*k+2-3*k)/2>=k]; print('P k',okP); assert okP==[2]
# Q: (CHO)_{2k}, k COOH, MW 58k in 100-150
okQ=[k for k in range(1,9) if (2*k)%2==0 and (2*2*k+2-2*k)/2>=k and 100<=58*k<=150]; print('Q k',okQ); assert okQ==[2]
# 物質の確認（原子量 H=1.0,C=12,O=16 で組成）
for smi,f,c,h in [('OC(=O)CCC(=O)O','C4H6O4',40.7,5.1),('OC(=O)/C=C\\C(=O)O','C4H4O4',41.4,3.4)]:
    m=Chem.MolFromSmiles(smi); assert CalcMolFormula(m)==f
    import re
    n={e:int(x or 1) for e,x in re.findall(r'([A-Z][a-z]?)(\d*)',f)}
    M=12*n['C']+1.0*n['H']+16*n['O']; print(f,M,round(12*n['C']/M*100,1),round(n['H']/M*100,1))
    assert round(12*n['C']/M*100,1)==c and round(n['H']/M*100,1)==h
# C4H6O4 の2価カルボン酸の異性体は2つ（コハク酸・メチルマロン酸）
cands=['OC(=O)CCC(=O)O','CC(C(=O)O)C(=O)O']
assert all(CalcMolFormula(Chem.MolFromSmiles(s))=='C4H6O4' for s in cands)
# 酸無水物の環の員数: コハク酸無水物は五員環
for s in ['O=C1CCC(=O)O1','CC1C(=O)OC1=O']:
    m=Chem.MolFromSmiles(s); print(s,[len(r) for r in m.GetRingInfo().AtomRings()])
print('OK')
