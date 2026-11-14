# ORG-STR-072 検証：元素分析・水素付加・イソプレン則・酸化開裂で候補 A〜E から X を決める
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import combinations
cand={'A':'CC1=CCC(CC1)C(=C)C','B':'C1CCC2=C(C1)CCCC2','C':'C=C1CCC2CC1C2(C)C','D':'CC1=CCC2CC1C2(C)C','E':'CC(=CCCC(=C)C=C)C'}
# (1) 元素分析
mC=22.0*12/44; mH=7.20*2/18
print('C,H mg',mC,mH,'合計',mC+mH)
assert abs(mC+mH-6.80)<1e-9
nC=mC/12; nH=mH/1.0
print('C:H',nC,nH,nH/nC)  # 0.5:0.8 = 5:8
assert abs(nH/nC-1.6)<1e-9
M=136.0; nX=6.80/M; nH2=1.12/22.4
print('X mmol',nX,'H2 mmol',nH2,'比',nH2/nX)
assert abs(nH2/nX-1)<1e-9
def cg(smi):
    m=Chem.MolFromSmiles(smi)
    C=[a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='C']
    E=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in m.GetBonds()]
    return m,C,E
def isopentane(S,E):
    es=[e for e in E if e[0] in S and e[1] in S]
    for T in combinations(es,4):
        deg={v:0 for v in S}; adj={v:[] for v in S}
        for a,b in T: deg[a]+=1;deg[b]+=1;adj[a].append(b);adj[b].append(a)
        seen={next(iter(S))};st=list(seen)
        while st:
            v=st.pop()
            for w in adj[v]:
                if w not in seen: seen.add(w);st.append(w)
        if len(seen)==5 and sorted(deg.values())==[1,1,1,2,3]:
            c3=[v for v in S if deg[v]==3][0]; c2=[v for v in S if deg[v]==2][0]
            if c2 in adj[c3]: return True
    return False
def isoprene(smi):
    m,C,E=cg(smi)
    return any(isopentane(set(S),E) and isopentane(set(C)-set(S),E) for S in combinations(C,5))
def cleave(smi):
    m=Chem.RWMol(Chem.MolFromSmiles(smi))
    dbl=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in m.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE]
    for a,b in dbl: m.RemoveBond(a,b)
    frags=Chem.GetMolFrags(m.GetMol(),asMols=False,sanitizeFrags=False)
    return [len(f) for f in frags]
res={}
for k,s in cand.items():
    m=Chem.MolFromSmiles(s)
    f=CalcMolFormula(m)
    ndb=sum(1 for b in m.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE)
    rings=m.GetRingInfo().NumRings()
    res[k]=(f,ndb,rings,isoprene(s),cleave(s))
    print(k,res[k])
    assert f=='C10H16'
ok=[k for k,v in res.items() if v[1]==1]
print('水素1分子:',ok); assert ok==['B','C','D']
ok=[k for k in ok if res[k][3]]
print('イソプレン則:',ok); assert ok==['C','D']
ok=[k for k in ok if res[k][4]==[10]]
print('開裂で1分子・炭素10:',ok); assert ok==['D']
# B の開裂生成物も炭素10の1分子（イソプレン則がないと区別できないこと）
assert res['B'][4]==[10] and res['C'][4]==[1,9]
# A, E はイソプレン則を満たす（水素の量で除かれる）
assert res['A'][3] and res['E'][3]
print('OK')
