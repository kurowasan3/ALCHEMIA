# ORG-ALC-046 検証：直鎖ヘプタンテトラオールの過ヨウ素酸酸化（消費量・生成物）の全異性体の総当たり
from itertools import combinations
from rdkit import Chem
def frag_smiles(n, S, cuts):
    # chain C1..Cn, OH at S, cut bonds (i,i+1) in cuts. Returns list of fragment SMILES after periodate oxidation:
    # each carbon at a cut end that had OH becomes C=O (aldehyde/ketone/acid if both sides cut -> HCOOH for CH(OH) with both bonds cut)
    segs=[]; start=1
    for i in range(1,n):
        if (i,i+1) in cuts: segs.append((start,i)); start=i+1
    segs.append((start,n))
    out=[]
    for a,b in segs:
        atoms=[]
        for c in range(a,b+1):
            cut_left = (c-1,c) in cuts; cut_right=(c,c+1) in cuts
            ncut=cut_left+cut_right
            atoms.append((c,ncut))
        # build SMILES linear
        s=''
        for k,(c,ncut) in enumerate(atoms):
            if ncut==0: sub='O' if c in S else ''
            elif ncut==1: sub='=O'
            else: sub='(=O)O'  # carbon cut on both sides -> HCOOH
            s+='C'+('('+sub+')' if sub and k<len(atoms)-1 else sub)
        m=Chem.MolFromSmiles(s); out.append(Chem.MolToSmiles(m))
    return sorted(out)
def periodate(n,S):
    cuts=[(i,i+1) for i in range(1,n) if i in S and i+1 in S]
    return len(cuts), frag_smiles(n,S,set(cuts))
def isomers(n,k):
    seen=set(); res=[]
    for S in combinations(range(1,n+1),k):
        R=tuple(sorted(n+1-i for i in S))
        key=min(S,R)
        if key in seen: continue
        seen.add(key); res.append(key)
    return res
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit.Chem import FindMolChiralCenters
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert check_equation('C7H16O4 + 2IO4- -> HCHO + HCOOH + C5H10O2 + 2IO3- + H2O') is None
# 数値：試料 0.328 g（C7H16O4 = 164）
M=7*12+16*1.0+4*16; n=0.328/M; print('M',M,'mmol',n*1000)
assert abs(n-2.00e-3)<1e-9
nIO4=4.00e-3; nHCOOH=0.100*20.0e-3
print('IO4/分子',nIO4/n,'HCOOH/分子',nHCOOH/n)
def iodoform(smi):
    m=Chem.MolFromSmiles(smi)
    return m.HasSubstructMatch(Chem.MolFromSmarts('[CH3][CX3](=O)')) or m.HasSubstructMatch(Chem.MolFromSmarts('[CH3][CX4H]([OH])'))
cands=[]
for S in isomers(7,4):
    k,fr=periodate(7,set(S))
    hcho=fr.count('C=O'); hcooh=fr.count('O=CO'); others=[f for f in fr if f not in ('C=O','O=CO')]
    if (k,hcho,hcooh)==(2,1,1): cands.append((S,others))
print('消費量・HCHO・HCOOH が合う候補',cands)
assert [c[0] for c in cands]==[(1,2,3,5),(1,2,3,6),(1,2,3,7)]
hit=[c for c in cands if any(iodoform(o) for o in c[1])]
print('ヨードホルム反応を示す生成物をもつもの',hit)
assert len(hit)==1 and hit[0][0]==(1,2,3,6) and hit[0][1]==['CC(O)CCC=O']
# 生成物 4-ヒドロキシペンタナールの分子式
assert CalcMolFormula(Chem.MolFromSmiles('CC(O)CCC=O'))=='C5H10O2'
# 1,2,3,6-ヘプタンテトラオールの不斉炭素と立体異性体
t=Chem.MolFromSmiles('OCC(O)C(O)CCC(O)C')
cc=FindMolChiralCenters(t,includeUnassigned=True,useLegacyImplementation=False); print('不斉炭素',cc)
assert len(cc)==3
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
st=set(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(t)); print('立体異性体',len(st)); assert len(st)==8
# 他の候補の生成物：3-ヒドロキシペンタナール・5-ヒドロキシペンタナールはヨードホルム陰性
assert not iodoform('CCC(O)CC=O') and not iodoform('O=CCCCCO')
print('OK')
