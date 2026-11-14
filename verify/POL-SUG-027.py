from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
units=R('3.24')/162; ag=R('0.432')/108; chains=ag/2
print(units, ag, chains, units/chains); assert units==R('0.02') and units/chains==10
assert units*2*108==R('4.32')
b=3-1; n=9; assert (b,n)==(2,9)
# build branched 9-mer: glucose ring template with attachment points
# glc(O4,O6 substituents): ring C1(O-anomeric) C2(OH) C3(OH) C4(O-x) C5 C6(O-y)
RING=[0]
def glc(o1, o4='O', o6='O'):
    RING[0]+=1; r=f"%{10+RING[0]}"
    return f"C{r}({o1})OC(C{o6})C({o4})C(O)C{r}O"
# build tree recursively: node = (child4, child6)
def build(node, root=False):
    c4, c6 = node
    o4 = 'O' + build(c4) if c4 is not None else 'O'
    o6 = 'O' + build(c6) if c6 is not None else 'O'
    o1 = 'O' if root else ''
    s = glc('X', o4, o6)
    return s.replace('(X)', '(O)') if root else s.replace('(X)', '', 1)
# chain of 5 with two 2-unit branches at units 2 and 4 -> 9 units, 2 alpha-1,6
leaf=(None,None)
br=((None,None),None)          # 2-unit branch: unit with a 1,4-linked child
tree=( ( ( ( (None,None), br), None), br), None)
def count(t): return 0 if t is None else 1+count(t[0])+count(t[1])
def c6(t): return 0 if t is None else (1 if t[1] is not None else 0)+c6(t[0])+c6(t[1])
print(count(tree), c6(tree)); assert count(tree)==9 and c6(tree)==2
smi=build(tree, root=True); print(smi); m=Chem.MolFromSmiles(smi)
assert Chem.rdMolDescriptors.CalcNumRings(m)==9 if hasattr(Chem,'rdMolDescriptors') else True
from rdkit.Chem import rdMolDescriptors as D; assert D.CalcNumRings(m)==9; f=CalcMolFormula(m); print(f); assert f=='C54H92O46'
print('POL-SUG-027 OK')
