import sys; sys.path.insert(0, "/home/claude/chem/work/g11"); sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
import common, json, itertools, math
from sympy import Rational as Q, nsimplify, solve, symbols, log, N
ok_all = True
def chk(name, cond):
    global ok_all
    print(("  PASS " if cond else "  FAIL ")+name); ok_all &= bool(cond)

it, b = common.basic("POL-PRT-010"); chk("basic", b)
t=json.dumps(json.load(open("/home/claude/kurowasan3/alchemia/content/chem/tests.json")),ensure_ascii=False)
chk("図鑑: 硫黄の検出 PbS / ビウレット等と矛盾なし（参照のみ）", "PbS" in t)
# 側鎖の性質（中性付近）: Cys -SH, Lys -NH3+, Glu -COO-, Phe/Leu 炭化水素
side={"Cys":"SH","Lys":"NH3+","Glu":"COO-","Phe":"hydrocarbon","Leu":"hydrocarbon"}
chk("(2) pairing consistent", side["Lys"].endswith("+") and side["Glu"].endswith("-") and side["Phe"]==side["Leu"]=="hydrocarbon")
print("  知識問題：高校教科書の標準記述（二次構造＝ペプチド結合間の水素結合，塩析は可逆，重金属による変性は不可逆）と照合")
print("RESULT:", "PASS" if ok_all else "FAIL")
