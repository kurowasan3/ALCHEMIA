import itertools
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
A={'Gly':75,'Ala':89,'Ser':105,'Cys':121,'Glu':147,'Lys':146,'Phe':165}
ac = 2*12+2*1+16; assert ac==42
M_X = R('2.93')*42/(R('3.35')-R('2.93')); print('M_X', M_X); assert M_X==293
sols=[c for c in itertools.combinations_with_replacement([a for a in A if a!='Lys'],3) if sum(A[x] for x in c)-36==293]
print(sols); assert sols==[('Gly','Ala','Phe')]
for a,exp in [(1,137),(2,274),(3,411)]:
    M=R('2.74')*42*a/(R('3.58')-R('2.74')); print(a,M); assert M==exp
assert 75*3-36==189 and 137<189
pairs=[c for c in itertools.combinations_with_replacement([a for a in A if a!='Lys'],2) if sum(A[x] for x in c)==274+36-146]
print(pairs); assert pairs==[('Gly','Ala')]
assert 411+36-2*146==155 and 155 not in A.values()
def nominal(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); w={'C':12,'H':1,'N':14,'O':16}
    return sum(w[a.GetSymbol()] for a in m.GetAtoms()), CalcMolFormula(m)
print(nominal('NCC(=O)NC(C)C(=O)NC(Cc1ccccc1)C(=O)O'))   # Gly-Ala-Phe
assert nominal('NCC(=O)NC(C)C(=O)NC(Cc1ccccc1)C(=O)O')[0]==293
print(nominal('NCCCCC(N)C(=O)NC(C)C(=O)NCC(=O)O'))       # Lys-Ala-Gly
assert nominal('NCCCCC(N)C(=O)NC(C)C(=O)NCC(=O)O')[0]==274
print('POL-PRT-027 OK')
