# ORG-CBX-017 検証：グリセリン・グリシン・酢酸からなる化合物の異性体
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g12')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item

GLY_E = 'C(=O)CN'          # グリシンのエステル（O-CO-CH2-NH2）
AC_E = 'C(=O)C'            # 酢酸エステル
NAG_E = 'C(=O)CNC(C)=O'    # N-アセチルグリシンのエステル
def glycerol(sub):  # sub: 位置1,2,3 の置換基（'' は OH の H）
    o = [('O' + s) if s else 'O' for s in sub]
    return f'C({o[0]})C({o[1]})C{o[2]}'
found = {}
for a, b in itertools.permutations(range(3), 2):      # 型 I
    sub = ['', '', '']; sub[a] = GLY_E; sub[b] = AC_E
    found.setdefault(Chem.CanonSmiles(glycerol(sub)), f'I gly@{a+1} ac@{b+1}')
for a in range(3):                                     # 型 II
    sub = ['', '', '']; sub[a] = NAG_E
    found.setdefault(Chem.CanonSmiles(glycerol(sub)), f'II @{a+1}')
chiral = 0
for smi, lab in found.items():
    m = Chem.MolFromSmiles(smi)
    cc = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    nh2 = m.HasSubstructMatch(Chem.MolFromSmarts('[NX3;H2]'))
    print(lab, smi, CalcMolFormula(m), '不斉', cc, 'NH2' if nh2 else '')
    assert CalcMolFormula(m) == 'C7H13NO5'
    chiral += bool(cc)
print('種類', len(found), '不斉あり', chiral)
assert len(found) == 5 and chiral == 4
Q = Chem.CanonSmiles('OCC(CO)OC(=O)CNC(C)=O')
assert Q in found and not Chem.FindMolChiralCenters(Chem.MolFromSmiles(Q), includeUnassigned=True, useLegacyImplementation=False)
assert (12*2 + 1.0*3 + 16) - 1.0 == 42 and 2 * 42 == 84
assert not check_item('ORG-CBX-017')
print('ORG-CBX-017 OK')
