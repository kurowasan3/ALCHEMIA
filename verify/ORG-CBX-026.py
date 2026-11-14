import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
dopc='CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC'
m=Chem.AddHs(Chem.MolFromSmiles(dopc)); f=CalcMolFormula(m); print(f); assert f=='C44H84NO8P'
w={'C':12,'H':1,'N':14,'O':16,'P':31}; ML=sum(w[a.GetSymbol()] for a in m.GetAtoms()); print(ML); assert ML==785==2*282+221
assert check_equation('C44H84NO8P + 4H2O -> C3H8O3 + 2C18H34O2 + H2PO4- + C5H14NO+') is None
M=symbols('M'); s=solve(R('7.85')*2/(2*M+221)-R('5.64')/M, M); print(s); assert s==[282]
sol=[(n,(14*n-250)//2) for n in range(10,26,2) if (14*n-250)>=0 and (14*n-250)%2==0]
print(sol); assert sol[0]==(18,1)
I2=R('0.0200')*254; print(float(I2)); assert I2==R('5.08')
assert R('7.85')/785==R('0.01')
print('ORG-CBX-026 OK')
