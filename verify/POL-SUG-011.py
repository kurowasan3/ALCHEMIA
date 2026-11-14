import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-SUG-011"); chk("basic", b)
t=json.load(open("/home/claude/kurowasan3/alchemia/content/chem/tests.json")); t=t if isinstance(t,list) else t.get("tests",t)
io=[x for x in t if isinstance(x,dict) and x.get("id")=="iodine-starch"][0]
chk("図鑑: 加熱で消え冷却で戻る", "加熱すると色が消え" in io["pos"]["starch"] and "戻る" in io["pos"]["starch"])
p=json.load(open("/home/claude/kurowasan3/alchemia/content/chem/polymers.json")); ps=json.dumps(p,ensure_ascii=False)
chk("図鑑: アミロース α-1,4 らせん / セルロース β-1,4", "らせん構造" in ps and "β-1,4-グリコシド結合" in ps and "α-1,4-グリコシド結合" in ps)
print("RESULT:", "PASS" if ok_all else "FAIL")
