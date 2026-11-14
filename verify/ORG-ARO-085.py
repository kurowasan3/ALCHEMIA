# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/ORG-ARO-085/item.json'
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

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ['src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer', 'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig']:
    assert k in it, k
# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
def MW(s): return Descriptors.MolWt(Chem.MolFromSmiles(s))
def F(s): return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s))
def mono_cl(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs() == 1:
            rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(17)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return len(out)
gallic = 'OC(=O)c1cc(O)c(O)c(O)c1'
D_m = 'OC(=O)c1cc(O)c(O)c(OC(=O)c2cc(O)c(O)c(O)c2)c1'
D_p = 'OC(=O)c1cc(O)c(OC(=O)c2cc(O)c(O)c(O)c2)c(O)c1'
for d in (D_m, D_p):
    assert F(d) == 'C14H10O9', F(d); print('D', round(MW(d), 1))
Dp_m = 'OC(=O)c1cc(OC)c(OC)c(OC(=O)c2cc(OC)c(OC)c(OC)c2)c1'
assert F(Dp_m) == 'C19H20O9'; print("D' MW", round(MW(Dp_m), 1))
Z = 'OC(=O)c1cc(OC)c(OC)c(OC)c1'; Y_m = 'OC(=O)c1cc(O)c(OC)c(OC)c1'; Y_p = 'OC(=O)c1cc(OC)c(O)c(OC)c1'
assert F(Z) == 'C10H12O5' and F(Y_m) == 'C9H10O5' and F(Y_p) == 'C9H10O5'
print('Z', round(MW(Z)), 'Y', round(MW(Y_m)))
assert round(MW(Z)) == 212 and round(MW(Y_m)) == 198
assert mono_cl(Y_m) == 2 and mono_cl(Y_p) == 1
# integer MW with H=1,C=12,O=16
def iw(c, h, o): return 12*c + h + 16*o
assert iw(7,6,5) == 170 and iw(14,10,9) == 322 and iw(19,20,9) == 392
for n in (2, 3, 4):
    print(n, round(170*n/(152*n+18), 4))
assert round(5.10/4.83, 4) == round(340/322, 4)
assert round(4.83/322*392, 2) == 5.88
print('OK')
