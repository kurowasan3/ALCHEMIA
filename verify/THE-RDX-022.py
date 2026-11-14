# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from fractions import Fraction as Fr
from sympy import symbols, solve
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
ITEM = '/home/claude/chem2/outh/items/THE-RDX-022/item.json'
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
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and ' + ' in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

def counts(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    c = {'C': 0, 'H': 0, 'O': 0}
    for a in m.GetAtoms(): c[a.GetSymbol()] += 1
    return c['C'], c['H'], c['O']
def mR(smi):
    x, y, z = counts(smi)
    m = Fr(2*z - y, x); R = Fr(2)*(4 - m)/3
    # 完全燃焼 CxHyOz + (x + y/4 - z/2) O2 -> x CO2 + y/2 H2O から直接
    nO2 = Fr(x) + Fr(y, 4) - Fr(z, 2)
    R2 = nO2*32/(12*x)
    assert R == R2, (smi, R, R2)
    return m, R
tab = {'メタノール': 'CO', 'ギ酸': 'OC=O', 'グルコース': 'OCC1OC(O)C(O)C(O)C1O', '酢酸': 'CC(=O)O',
       'グリコール酸': 'OCC(=O)O', 'グリオキシル酸': 'O=CC(=O)O', 'シュウ酸': 'OC(=O)C(=O)O', 'メタン': 'C', 'エタン': 'CC'}
res = {k: mR(v) for k, v in tab.items()}
for k, v in res.items(): print(k, CalcMolFormula(Chem.MolFromSmiles(tab[k])), 'm =', v[0], 'R =', v[1])
assert res['メタノール'][1] == 4 and res['ギ酸'][1] == Fr(4, 3) and res['グルコース'][1] == Fr(8, 3) and res['メタン'][1] == Fr(16, 3)
assert max(r for _, r in res.values()) == Fr(16, 3)
# (3)
R = Fr(128, 48); m = 4 - Fr(3, 2)*R; print('(3) m', m); assert m == 0
nC = Fr(48, 12)/1000; nacid = Fr(2, 1000); print('C/分子', nC/nacid); assert nC/nacid == 2
cands = ['CC(=O)O', 'OCC(=O)O', 'O=CC(=O)O']   # C2 で COOH 1 個
ok = [s for s in cands if mR(s)[0] == 0]; print('m=0 の候補', ok); assert ok == ['CC(=O)O']
assert 2*Fr(2, 1000)*32*1000 == 128   # 酢酸 2.0 mmol は O2 4.0 mmol = 128 mg
# (4)
a, b = symbols('a b')
s = solve([a + b - 3.0e-3, 6*a + 2*b - 4*128e-3/32], [a, b]); print('(4)', s, s[a]*32e3, s[b]*46e3)
assert round(float(s[a]*32e3)) == 80 and round(float(s[b]*46e3)) == 23
# (5)
mbar = (res['メタノール'][0] + res['ギ酸'][0])/2; print('(5) 1:1 平均 m', mbar); assert mbar == 0 and Fr(2)*(4 - mbar)/3 == Fr(8, 3)
print('全検証OK:', ITEM)
