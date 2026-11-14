# POL-PRT-072 ジニトロフルオロベンゼンによるアミノ基の標識
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import permutations
W={'C':12,'H':1,'N':14,'O':16,'F':19}
def mw(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); return sum(W[a.GetSymbol()] for a in m.GetAtoms())
DNP='c1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]'
mols={'FDNB':'Fc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]',
 'DNP-Gly':'OC(=O)CNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]',
 'bisDNP-Lys':'OC(=O)C(Nc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-])CCCCNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]',
 'eDNP-Lys':'OC(=O)C(N)CCCCNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]'}
for k,s in mols.items(): print(k,CalcMolFormula(Chem.MolFromSmiles(s)),mw(s))
# 配列の場合分け：Lys(K) Ala(A) Gly(G) 各1
def labeled(seq):
    out=[]
    for i,a in enumerate(seq):
        n=(1 if i==0 else 0)+(1 if a=='K' else 0)
        out.append((a,n))
    return out
for seq in permutations('KAG'):
    lab=labeled(seq); print(''.join(seq),lab,'DNP total',sum(n for _,n in lab))
# 条件X: 2個標識の化合物1つ+未標識2 → Kが N末端
print('X:',[''.join(s) for s in permutations('KAG') if sorted(n for _,n in labeled(s))==[0,0,2]])
print('Y:',[''.join(s) for s in permutations('KAG') if dict(labeled(s))=={'A':1,'K':1,'G':0}])
