from rdkit import Chem
POS={'T':(0,1),'R':(0.95,-0.1),'L':(-0.95,-0.1),'B':(0,-1)}
grp={'COOH':('C',[('O',2),('O',1)]),'NH2':('N',[]),'H':('H',[]),'CH3':('C',[]),'CH2SH':('C',[('S',1)])}
def build(assign):
    atoms=[('C',0,0)]; bonds=[]
    for p,g in assign.items():
        x,y=POS[p]; sym,subs=grp[g]; idx=len(atoms); atoms.append((sym,x,y))
        bonds.append((1,idx+1,1,{'L':6,'B':1}.get(p,0)))
        for k,(s2,o) in enumerate(subs):
            j=len(atoms); atoms.append((s2,x*1.8+0.3*k,y*1.8+0.3*k)); bonds.append((idx+1,j+1,o,0))
    mb="\n  RDKit          2D\n\n%3d%3d  0  0  0  0  0  0  0  0999 V2000\n"%(len(atoms),len(bonds))
    for s,x,y in atoms: mb+="%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n"%(x,y,0,s)
    for a,b,o,f in bonds: mb+="%3d%3d%3d%3d\n"%(a,b,o,f)
    mb+="M  END\n"
    m=Chem.MolFromMolBlock(mb,removeHs=False)
    return Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(m)))
Lala=Chem.CanonSmiles('C[C@@H](C(=O)O)N'); Dala=Chem.CanonSmiles('C[C@H](C(=O)O)N'); Lcys=Chem.CanonSmiles('C([C@@H](C(=O)O)N)S')
order=['T','R','L','B']
def parity(base, new):
    perm=[base.index(g) for g in new]; s=0
    for i in range(4):
        for j in range(i+1,4):
            if perm[i]>perm[j]: s+=1
    return s%2
base=['COOH','NH2','H','CH3']
assert build(dict(zip(order,base)))==Lala
cases={'ア':(['NH2','COOH','H','CH3'],'D'),'イ':(['CH3','NH2','H','COOH'],'D'),'ウ':(['NH2','COOH','CH3','H'],'L'),'エ':(['NH2','CH3','H','COOH'],'L')}
for k,(arr,exp) in cases.items():
    smi=build(dict(zip(order,arr))); got='L' if smi==Lala else ('D' if smi==Dala else '?')
    par='D' if parity(base,arr) else 'L'
    print(k,smi,got,par); assert got==exp==par
o=['CH2SH','COOH','H','NH2']; smi=build(dict(zip(order,o))); print('オ',smi); assert smi==Lcys
assert parity(['COOH','NH2','H','CH2SH'],o)==0
print('POL-PRT-031 OK')
