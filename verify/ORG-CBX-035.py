"""ORG-CBX-035 検証：C8H8O2 芳香族エステルの候補・NaOH 量・クレゾールの一塩素置換体の数・反応式"""
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Descriptors
esters = {'安息香酸メチル': 'COC(=O)c1ccccc1', '酢酸フェニル': 'CC(=O)Oc1ccccc1', 'ギ酸ベンジル': 'O=COCc1ccccc1',
          'ギ酸o-メチルフェニル': 'O=COc1ccccc1C', 'ギ酸m-メチルフェニル': 'O=COc1cccc(C)c1', 'ギ酸p-メチルフェニル': 'O=COc1ccc(C)cc1'}
for k, s in esters.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C8H8O2'; 
# 他に C8H8O2・ベンゼン環1個・エステル(-C(=O)O-C) で、環以外の炭素が2個のものがないか：置換基の組み合わせを列挙
smis = ['COC(=O)c1ccccc1', 'CC(=O)Oc1ccccc1', 'O=COCc1ccccc1', 'O=COc1ccccc1C', 'O=COc1cccc(C)c1', 'O=COc1ccc(C)cc1',
        'O=C(O)Cc1ccccc1']
patt = Chem.MolFromSmarts('[#6][CX3](=O)[OX2][#6]'); patt2 = Chem.MolFromSmarts('[CX3H1](=O)[OX2][#6]')
nest = 0
for s in smis:
    m = Chem.MolFromSmiles(s)
    if m.HasSubstructMatch(patt) or m.HasSubstructMatch(patt2): nest += 1
assert nest == 6

# 置換基を列挙して C8H8O2・ベンゼン環1個のエステルを数え上げる（ラクトン環は C8H8O2 で環1個の条件では生じない）
groups = ['C','CC','O','OC','C=O','CO','C(=O)O','OC=O','C(=O)OC','OC(C)=O','COC=O','CC(=O)O','OCC=O','CC(=O)OC','C(C)OC=O','COC(C)=O',
          'C(=O)C','OC(=O)C','C(=O)CO','CC=O','OCC','COC','CCO','C(O)C=O','C(=O)OCC','CCOC=O']
found = set()
def benz(subs):
    # subs: list of (position, group) positions 0..5
    ring = ['c'] * 6; tail = [''] * 6
    for pos, g in subs: tail[pos] = '(' + g + ')'
    smi = ''.join(('c' + tail[i]) if i else ('c1' + tail[0]) for i in range(6)) + '1'
    return smi
for g1 in groups:
    smi = benz([(0, g1)]); m = Chem.MolFromSmiles(smi)
    if m and CalcMolFormula(m) == 'C8H8O2' and (m.HasSubstructMatch(patt) or m.HasSubstructMatch(patt2)): found.add(Chem.MolToSmiles(m))
    for g2 in groups:
        for pos in (1, 2, 3):
            m = Chem.MolFromSmiles(benz([(0, g1), (pos, g2)]))
            if m and CalcMolFormula(m) == 'C8H8O2' and (m.HasSubstructMatch(patt) or m.HasSubstructMatch(patt2)): found.add(Chem.MolToSmiles(m))
print('列挙したエステル', len(found), sorted(found))
assert found == {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in esters.values()}

M = 8*12 + 8*1.0 + 2*16; assert abs(M-136) < 1e-9
n = 1.36/M; print('mol', n)
for g, exp in [(0.400, 1), (0.800, 2)]:
    assert abs((g/40.0)/n - exp) < 1e-9
# クレゾールの環の一塩素置換体の数
def ring_mono(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs() == 1:
            em = Chem.RWMol(m); i = em.AddAtom(Chem.Atom(17)); em.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
            em.GetAtomWithIdx(a.GetIdx()).SetNumExplicitHs(0)
            mm = em.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return len(out)
cnt = {k: ring_mono(s) for k, s in {'o': 'Cc1ccccc1O', 'm': 'Cc1cccc(O)c1', 'p': 'Cc1ccc(O)cc1'}.items()}
print('一塩素置換体', cnt); assert cnt == {'o': 4, 'm': 4, 'p': 2}
it = json.load(open('/home/claude/chem5/out/items/ORG-CBX-035/item.json'))
k = 0
for key in ('problem', 'answer', 'solution'):
    for m in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[key]):
        if re.search(r'->|<=>', m):
            r = check_equation(m); print(m, r or 'OK'); assert r is None; k += 1
assert k >= 4
print('ALL OK')
