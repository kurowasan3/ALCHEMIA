# 検証スクリプト（グループ g11）  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-RDX-055/item.json'
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
assert check_equation('2XeF2 + 2H2O -> 2Xe + 3HF + O2') is not None  # 検査器が誤りを検出できること
def sig(x, n):
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
import sympy as sp
a, b, c = sp.symbols('a b c')
Xe = a + sp.Rational(2, 3)*b; O2 = (a + b)/2; XeO3 = b/3 + c; HF = 2*a + 4*b + 6*c
assert sp.expand(4*O2 + 6*XeO3 - HF) == 0
def rank(exprs):
    return sp.Matrix([[sp.diff(e, v) for v in (a, b, c)] for e in exprs]).rank()
assert rank([O2, XeO3, HF]) == 2
mass = 169*a + 207*b + 245*c
F_amount = HF
print('ア', rank([O2, XeO3, HF, Xe]), 'イ', rank([O2, XeO3, HF, mass]), 'ウ', rank([O2, XeO3, HF, F_amount]))
assert rank([O2, XeO3, HF, Xe]) == 3 and rank([O2, XeO3, HF, mass]) == 3 and rank([O2, XeO3, HF, F_amount]) == 2
# 試料 X
ngas = 59.4/22.4; nO2 = ngas*0.396; nXe = ngas - nO2; nXeO3 = 0.100*84.0/6
s = sp.solve([sp.Eq(O2, nO2), sp.Eq(Xe, nXe), sp.Eq(XeO3, nXeO3)], [a, b, c])
print('X', s); tot = sum(s.values())
pct = [round(float(s[v]/tot*100)) for v in (a, b, c)]; print(pct); assert pct == [20, 50, 30]
# 試料 Y
s2 = sp.solve([sp.Eq(O2, 0.800), sp.Eq(XeO3, 1.00), sp.Eq(mass, 463.0)], [a, b, c]); print('Y', s2)
assert abs(s2[a]-0.40) < 0.005 and abs(s2[b]-1.20) < 0.005 and abs(s2[c]-0.60) < 0.005
assert abs(float(HF.subs(s2)) - 9.20) < 1e-6
s3 = sp.solve([sp.Eq(O2, 0.800), sp.Eq(XeO3, 1.00), sp.Eq(mass, 464.0)], [a, b, c]); print('質量+1mg', s3[b]-s2[b])
assert abs(abs(s3[b]-s2[b]) - 0.0229) < 0.001
# 縮退方向 (a,b,c)->(a-3t,b+3t,c-t)
t = sp.symbols('t'); sub = {a: a-3*t, b: b+3*t, c: c-t}
assert sp.expand(O2.subs(sub, simultaneous=True) - O2) == 0 and sp.expand(XeO3.subs(sub, simultaneous=True) - XeO3) == 0
assert sp.expand(207 - 169 - sp.Rational(245, 3)) == sp.Rational(-131, 3)
print('OK')
