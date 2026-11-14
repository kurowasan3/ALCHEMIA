# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/outh/items/ORG-CBX-038/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from fractions import Fraction
AM = dict(H=1.0, C=12.0, N=14.0, O=16.0, Na=23.0, K=39.0, Cu=63.5, I=127.0)
from chemlib import parse_formula
def M(f):
    c, _ = parse_formula(f)
    return sum(AM[e] * float(n) for e, n in c.items())

import sympy as sp
TG, DG, MG = M('C51H98O6'), M('C35H68O5'), M('C19H38O4')
assert (TG, DG, MG) == (806, 568, 330) and M('C3H8O3') == 92
assert TG - DG == 238 and DG - MG == 238 and M('C16H31O') - M('H') == 238
t, d, m, g, m2 = sp.symbols('t d m g m2')
sol = sp.solve([t+d+m+g-sp.Rational(1,10), d+2*m+3*g-sp.Rational(17,100), 2*g+m-sp.Rational(95,1000), g-sp.Rational(4,100)], [t,d,m,g])
assert (sol[t], sol[d], sol[m]) == (sp.Rational(25,1000), sp.Rational(20,1000), sp.Rational(15,1000))
# (c) は (a) と恒等的に同じ，(d) は (a),(b) から決まる
mass = sol[t]*TG + sol[d]*DG + sol[m]*MG
assert abs(float(mass) - (806*0.100 - 238*0.170 - 92*0.0400)) < 1e-9 and sig(float(mass),3) == 36.5
tt, dd, mm = sp.symbols('tt dd mm')
expr = (tt*TG + dd*DG + mm*MG) - (TG*(tt+dd+mm) - 238*(dd+2*mm))
assert sp.simplify(expr) == 0
# 2-MG を含む場合
m1 = sp.Rational(15,1000)
s2 = sp.solve([t+d+m1+m2+g-sp.Rational(1,10), d+2*(m1+m2)+3*g-sp.Rational(17,100), g-sp.Rational(4,100)], [t,d,g])
assert sp.simplify(s2[t] - (sp.Rational(25,1000)+m2)) == 0 and sp.simplify(s2[d] - (sp.Rational(20,1000)-2*m2)) == 0
print('TG の範囲', s2[t].subs(m2,0), '〜', s2[t].subs(m2, sp.Rational(1,100)))
# 隣接ジオールの有無（RDKit）
from rdkit import Chem
patt = Chem.MolFromSmarts('[CX4;!$(C=O)]([OX2H])[CX4]([OX2H])')
R = 'CCCCCCCCCCCCCCCC(=O)'
mols = {'グリセリン': ('OCC(O)CO', 2), '1-MG': (R+'OCC(O)CO', 1), '2-MG': ('OCC(OC(=O)CCCCCCCCCCCCCCC)CO', 0),
        '1,2-DG': (R+'OCC(OC(=O)CCCCCCCCCCCCCCC)CO', 0), '1,3-DG': (R+'OCC(O)COC(=O)CCCCCCCCCCCCCCC', 0),
        'TG': (R+'OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC', 0)}
for k, (smi, want) in mols.items():
    mol = Chem.MolFromSmiles(smi)
    # 切れる C–C 結合の数＝OH をもつ炭素どうしの結合の数
    cuts = set()
    for a, _, b, _ in [mt for mt in mol.GetSubstructMatches(patt)]:
        cuts.add(tuple(sorted((a, b))))
    assert len(cuts) == want, (k, len(cuts))
print('ORG-CBX-038 数値・構造OK')
