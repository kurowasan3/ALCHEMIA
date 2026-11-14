import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
it, b = common.basic("POL-PRT-012"); chk("basic", b)
from rdkit.Chem import AllChem
AAS={"D-Ala":"C[C@@H](N)C(=O)O","L-Ala":"C[C@H](N)C(=O)O","L-Phe":"N[C@@H](Cc1ccccc1)C(=O)O"}
R=AAS
rxn=AllChem.ReactionFromSmarts("[C:1](=[O:2])[OX2H1].[NX3H2:3]>>[C:1](=[O:2])[N:3]")
def dip(a,b):
    ps={Chem.MolToSmiles(p[0]) for p in rxn.RunReactants((Chem.MolFromSmiles(AAS[a]),Chem.MolFromSmiles(AAS[b])))}
    good=[x for x in ps if Chem.MolFromSmiles(x).HasSubstructMatch(Chem.MolFromSmarts("[NX3H2]C(*)C(=O)NC(*)C(=O)[OH]"))]
    assert len(set(good))==1, (a,b,ps); return good[0]
def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace("@@","!").replace("@","@@").replace("!","@")))
chk("D/L-Ala are mirror images", mirror(AAS["L-Ala"])==Chem.MolToSmiles(Chem.MolFromSmiles(AAS["D-Ala"])))
from rdkit.Chem import FindMolChiralCenters
chk("L-Ala = S, L-Phe = S (CIP)", FindMolChiralCenters(Chem.MolFromSmiles(AAS["L-Ala"]))[0][1]=="S" and FindMolChiralCenters(Chem.MolFromSmiles(AAS["L-Phe"]))[0][1]=="S")
names=list(R); D={(a,b):dip(a,b) for a in names for b in names}
chk("(1) 9 distinct", len(set(D.values()))==9)
S=set(D.values()); pairs={frozenset((s,mirror(s))) for s in S if mirror(s) in S and mirror(s)!=s}
print(" pairs:",len(pairs)); chk("(2) 2 pairs", len(pairs)==2)
chk("no achiral dipeptide", all(mirror(s)!=s for s in S))
rac=[k for k in D if sorted(k)==["D-Ala","L-Ala"]]; print(" racemic hydrolysate:",rac)
chk("(3) 2 kinds racemic", len(rac)==2)
chk("(4) D-Ala-L-Ala optically active (≠ mirror)", mirror(D[("D-Ala","L-Ala")])!=D[("D-Ala","L-Ala")] and mirror(D[("D-Ala","L-Ala")])==D[("L-Ala","D-Ala")])
print("RESULT:", "PASS" if ok_all else "FAIL")
