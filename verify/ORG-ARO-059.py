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
d = check_all_equations('/home/claude/chem5/out/items/ORG-ARO-059/item.json')
# R = C7H7 でベンゼン環を含む基を総当たり（R-H = C7H8 の H を1つ外す）
toluene=Chem.MolFromSmiles('Cc1ccccc1')
Rs=set()
for a in toluene.GetAtoms():
    if a.GetTotalNumHs()==0: continue
    rw=Chem.RWMol(toluene); i=rw.AddAtom(Chem.Atom(0)); rw.AddBond(a.GetIdx(),i,Chem.BondType.SINGLE)
    m=rw.GetMol(); Chem.SanitizeMol(m); Rs.add(Chem.MolToSmiles(m))
print('R groups',Rs); assert len(Rs)==4   # o,m,p-tolyl + benzyl
Xs=set()
for r in Rs:
    r1=Chem.MolFromSmiles(r.replace('*','[*:1]')); r2=Chem.MolFromSmiles(r.replace('*','[*:2]'))
    core=Chem.MolFromSmiles('O=C(OC[*:2])[*:1]')
    combo=Chem.molzip(Chem.molzip(core,r1),r2)
    smi=Chem.MolToSmiles(combo); Xs.add(smi)
    assert CalcMolFormula(combo)=='C16H16O2', smi
    q=Chem.molzip(Chem.MolFromSmiles('OC[*:1]'),r1)   # Q = R-CH2OH
    assert q.HasSubstructMatch(Chem.MolFromSmarts('[CH2][OH]'))
print(len(Xs)); [print(' ',x) for x in sorted(Xs)]
assert len(Xs)==4
# 第二級アルコール（1-フェニルエタノール）を使うものは Q→ケトンなので除外される: 確認として C16H16O2 の式にもなる
assert CalcMolFormula(Chem.MolFromSmiles('CC(OC(=O)Cc1ccccc1)c1ccccc1'))=='C16H16O2'
# S の候補と環 H の種類
def ringH_classes(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); rk=list(Chem.CanonicalRankAtoms(m,breakTies=False))
    return len({rk[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol()=='H' and a.GetNeighbors()[0].GetIsAromatic()})
for name,s in [('phthalic','OC(=O)c1ccccc1C(=O)O'),('isophthalic','OC(=O)c1cccc(C(=O)O)c1'),('terephthalic','OC(=O)c1ccc(C(=O)O)cc1')]:
    print(name,ringH_classes(s))
assert ringH_classes('OC(=O)c1ccc(C(=O)O)cc1')==1 and ringH_classes('OC(=O)c1cccc(C(=O)O)c1')==3 and ringH_classes('OC(=O)c1ccccc1C(=O)O')==2
print('OK')
