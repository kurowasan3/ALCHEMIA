import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

from rdkit import Chem
it, b = common.basic("POL-SUG-017"); chk("basic", b)
def periodate(smi):
    m=Chem.MolFromSmiles(smi)
    def hasOH(a): return a.GetSymbol()=="C" and any(nb.GetSymbol()=="O" and nb.GetTotalNumHs()==1 for nb in a.GetNeighbors())
    cut=[bd for bd in m.GetBonds() if bd.GetBeginAtom().GetSymbol()=="C" and bd.GetEndAtom().GetSymbol()=="C" and hasOH(bd.GetBeginAtom()) and hasOH(bd.GetEndAtom())]
    deg={}
    for bd in cut:
        for a in (bd.GetBeginAtomIdx(),bd.GetEndAtomIdx()): deg[a]=deg.get(a,0)+1
    return len(cut), sum(1 for v in deg.values() if v==2)
cases={"グリセリン":"OCC(O)CO","エチレングリコール":"OCCO","1,3-プロパンジオール":"OCCCO",
 "非還元末端(4位OH)":"OCC1OC(OC)C(O)C(O)C1O","鎖の途中(4-O-置換)":"OCC1OC(OC)C(O)C(O)C1OC","分岐(4,6-O-置換)":"COCC1OC(OC)C(O)C(O)C1OC",
 "α-メチルグルコシド":"OCC1OC(OC)C(O)C(O)C1O"}
r={k:periodate(v) for k,v in cases.items()}
for k,v in r.items(): print("  ",k,"HIO4,HCOOH =",v)
chk("glycerol 2/1, EG 1/0, 1,3-PD 0/0", r["グリセリン"]==(2,1) and r["エチレングリコール"]==(1,0) and r["1,3-プロパンジオール"]==(0,0))
chk("only non-reducing end gives HCOOH", r["非還元末端(4位OH)"][1]==1 and r["鎖の途中(4-O-置換)"][1]==0 and r["分岐(4,6-O-置換)"][1]==0)
mol=Q(162,100)/(Q(324,100)*10**5); ends=Q(4,10)*Q(1,10**3)/mol; print(mol,ends)
chk("(3) 80 ends, 79 branches", ends==80 and ends-1==79)
print("RESULT:", "PASS" if ok_all else "FAIL")
