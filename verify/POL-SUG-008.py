import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-SUG-008"); chk("basic", b)
P = {"グルコース":dict(I=0,red=1,hyd=None),"フルクトース":dict(I=0,red=1,hyd=None),
     "スクロース":dict(I=0,red=0,hyd={"グルコース","フルクトース"}),"マルトース":dict(I=0,red=1,hyd={"グルコース"}),
     "ラクトース":dict(I=0,red=1,hyd={"グルコース","ガラクトース"}),"デンプン":dict(I=1,red=0,hyd={"グルコース"})}
sols=[]
for perm in itertools.permutations(P):
    m=dict(zip("ABCDEF",perm))
    if [k for k in "ABCDEF" if P[m[k]]["I"]]!=["F"]: continue
    if set(k for k in "ABCDEF" if P[m[k]]["red"])!=set("ABDE"): continue
    if P[m["D"]]["hyd"]!={m["A"]} or P[m["E"]]["hyd"]!={m["A"],"ガラクトース"}: continue
    sols.append(m)
print(sols); chk("unique assignment", len(sols)==1 and sols[0]["C"]=="スクロース" and sols[0]["B"]=="フルクトース")
s=sols[0]; chk("(2) only C yields fructose", [k for k in "ABCDEF" if P[s[k]]["hyd"] and s["B"] in P[s[k]]["hyd"]]==["C"])
t=json.load(open("/home/claude/kurowasan3/alchemia/content/chem/tests.json"))
t=t if isinstance(t,list) else t.get("tests",t)
f=[x for x in t if isinstance(x,dict) and x.get("id")=="fehling"][0]
chk("図鑑 fehling: glucose,fructose,maltose pos; sucrose,starch not pos", all(k in f["pos"] for k in ["glucose","fructose","maltose"]) and "sucrose" not in f["pos"] and "starch" not in f["pos"])
print("RESULT:", "PASS" if ok_all else "FAIL")
