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
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog('rdApp.*')
d = check_all_equations('/home/claude/chem5/outh/items/ORG-STR-046/item.json')
# --- (1) 総当たり：ベンゼンに C・O を1原子ずつ（単結合）または =O を付け足して木を成長させる
def grow(smi):
    m=Chem.MolFromSmiles(smi); out=set()
    for a in m.GetAtoms():
        if a.GetTotalNumHs()==0: continue
        for el,bt in ((6,Chem.BondType.SINGLE),(8,Chem.BondType.SINGLE),(8,Chem.BondType.DOUBLE)):
            if bt==Chem.BondType.DOUBLE and (a.GetIsAromatic() or a.GetSymbol()!='C' or a.GetTotalNumHs()<2): continue
            rw=Chem.RWMol(m); i=rw.AddAtom(Chem.Atom(el)); rw.AddBond(a.GetIdx(),i,bt)
            try:
                mm=rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
            except Exception: pass
    return out
lvl={'c1ccccc1'}; allm=set(lvl)
for _ in range(5):
    new=set()
    for s in lvl: new|=grow(s)
    lvl=new; allm|=new
ester=Chem.MolFromSmarts('[CX3;!R](=O)[OX2][#6]'); oo=Chem.MolFromSmarts('[#8]~[#8]')
E=sorted(s for s in allm if CalcMolFormula(Chem.MolFromSmiles(s))=='C9H10O2'
         and Chem.MolFromSmiles(s).HasSubstructMatch(ester) and not Chem.MolFromSmiles(s).HasSubstructMatch(oo))
print('total',len(E)); assert len(E)==24
def acid_class(s):
    m=Chem.MolFromSmiles(s)
    if m.HasSubstructMatch(Chem.MolFromSmarts('[CH1](=O)O')): return 'formate'
    if m.HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)O')): return 'acetate'
    if m.HasSubstructMatch(Chem.MolFromSmarts('[CH3][CH2]C(=O)O')): return 'propanoate'
    return 'aromatic acid'
from collections import Counter
cnt=Counter(acid_class(s) for s in E); print(cnt)
assert cnt=={'formate':14,'acetate':4,'propanoate':1,'aromatic acid':5}
# --- (2) フェノール類のエステル（O がベンゼン環に直接結合）
phen=Chem.MolFromSmarts('C(=O)Oc'); P=[s for s in E if Chem.MolFromSmiles(s).HasSubstructMatch(phen)]
print('phenyl esters',len(P)); assert len(P)==13
# --- (3) 不斉炭素
ch=[s for s in E if Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True,useLegacyImplementation=False)]
print('chiral',ch); assert len(ch)==1
# --- (4) 絞り込み
a_b=[s for s in P if acid_class(s)=='formate']; assert len(a_b)==9
def ringH_classes(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); rk=list(Chem.CanonicalRankAtoms(m,breakTies=False))
    return len({rk[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol()=='H' and a.GetNeighbors()[0].GetIsAromatic()})
c=[s for s in a_b if ringH_classes(s)==2]; print('after c',c); assert len(c)==3
def free_op(smi):   # フェノール Y の OH のオルト・パラ位で H のある位置の数
    rxn=Chem.MolFromSmiles(smi)
    m=Chem.MolFromSmiles(Chem.MolToSmiles(rxn).replace('OC=O','O'))
    o=[a for a in m.GetAtoms() if a.GetSymbol()=='O'][0]; c1=o.GetNeighbors()[0]
    dm=Chem.GetDistanceMatrix(m); n=0
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs()>0 and int(dm[c1.GetIdx()][a.GetIdx()]) in (1,3): n+=1
    return n
res={s:free_op(s) for s in c}; print(res)
X=[s for s in c if res[s]==1]; print('X',X)
assert X==[Chem.MolToSmiles(Chem.MolFromSmiles('O=COc1c(C)cccc1C'))]
print('OK')
