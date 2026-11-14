# POL-RES-013 検証：C5H8O2 のエステルの特定
import sys, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from g12lib import check_item
# 酸部分 R1-CO-（SMILES の先頭に C(=O) が続く形）とアルコール部分 -O-R2
ACIDS = {'ギ酸': ('[H]', 'O=CO'), '酢酸': ('C', 'CC(=O)O'), 'プロピオン酸': ('CC', 'CCC(=O)O'), '酪酸': ('CCC', 'CCCC(=O)O'),
         'イソ酪酸': ('CC(C)', 'CC(C)C(=O)O'), 'アクリル酸': ('C=C', 'C=CC(=O)O'), 'クロトン酸': ('CC=C', 'CC=CC(=O)O'),
         'メタクリル酸': ('C=C(C)', 'C=C(C)C(=O)O'), '3-ブテン酸': ('C=CC', 'C=CCC(=O)O')}
ALCS = {'メタノール': ('C', 'CO'), 'エタノール': ('CC', 'CCO'), '1-プロパノール': ('CCC', 'CCCO'), '2-プロパノール': ('C(C)C', 'CC(C)O'),
        'ビニル': ('C=C', None), 'アリル': ('CC=C', 'C=CCO'), '1-プロペニル': ('C=CC', None), 'イソプロペニル': ('C(=C)C', None),
        '1-ブタノール': ('CCCC', 'CCCCO'), '2-ブタノール': ('C(C)CC', 'CCC(C)O'), '2-メチル-1-プロパノール': ('CC(C)C', 'CC(C)CO'),
        'tert-ブタノール': ('C(C)(C)C', 'CC(C)(C)O'), '3-ブテニル': ('CCC=C', 'C=CCCO'), '2-ブテニル': ('CC=CC', 'CC=CCO'),
        '1-メチルアリル': ('C(C)C=C', 'C=CC(C)O'), 'メタリル': ('CC(C)=C', 'C=C(C)CO')}
def ctdb(m):
    return len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)}) > 1 and any(
        b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() != Chem.BondStereo.STEREONONE
        for x in EnumerateStereoisomers(m) for b in x.GetBonds())
iodo = Chem.MolFromSmarts('[CH3][CX4;H1,H2;!$(C(C)(C)C)]([OX2H1])')  # CH3CH(OH)-R または エタノール
iodo2 = Chem.MolFromSmarts('[CH3][CH2][OH]')
cands = []
for an, (r1, acid) in ACIDS.items():
    for ln, (r2, alc) in ALCS.items():
        smi = (f'{r1}C(=O)O{r2}' if r1 != '[H]' else f'O=CO{r2}')
        m = Chem.MolFromSmiles(smi)
        if CalcMolFormula(m) != 'C5H8O2': continue
        ncc = sum(b.GetBondType() == Chem.BondType.DOUBLE and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C' for b in m.GetBonds())
        if ncc != 1: continue
        acid_unsat = '=' in r1
        cands.append(dict(smi=Chem.CanonSmiles(smi), acid=an, alc=ln, acid_unsat=acid_unsat, cistrans=ctdb(Chem.MolFromSmiles(smi)),
                          alc_smi=alc, acid_C=Chem.MolFromSmiles(acid).GetNumAtoms() - 2 if acid else None,
                          branched=any(a.GetSymbol() == 'C' and sum(n.GetSymbol() == 'C' for n in a.GetNeighbors()) >= 3 for a in m.GetAtoms())))
cands = [c for c in cands if not c['cistrans']]
for c in cands: print(c['smi'], c['acid'], c['alc'], '酸側二重結合' if c['acid_unsat'] else 'アルコール側二重結合', '枝分かれ' if c['branched'] else '')
def iodoform(alc):
    m = Chem.MolFromSmiles(alc); return m.HasSubstructMatch(iodo) or m.HasSubstructMatch(iodo2)
P = [c for c in cands if c['acid_unsat'] and c['alc_smi'] and iodoform(c['alc_smi'])]
Q = [c for c in cands if not c['acid_unsat'] and c['acid'] == 'プロピオン酸']
R = [c for c in cands if c['acid_unsat'] and c['alc'] == 'メタノール' and c['branched']]
print('P', [c['smi'] for c in P], 'Q', [c['smi'] for c in Q], 'R', [c['smi'] for c in R])
assert [c['smi'] for c in P] == [Chem.CanonSmiles('C=CC(=O)OCC')]
assert [c['smi'] for c in Q] == [Chem.CanonSmiles('C=COC(=O)CC')]
assert [c['smi'] for c in R] == [Chem.CanonSmiles('C=C(C)C(=O)OC')]
T = Chem.MolFromSmiles('O=C1CCCCO1'); assert CalcMolFormula(T) == 'C5H8O2'
ri = T.GetRingInfo(); assert ri.NumRings() == 1 and len(ri.AtomRings()[0]) == 6
assert sum(T.GetAtomWithIdx(i).GetSymbol() == 'O' for i in ri.AtomRings()[0]) == 1
F = lambda s: re.sub(r'\*\d*', '', CalcMolFormula(Chem.MolFromSmiles(s)))
assert F('*CC(*)C(=O)O') == 'C3H4O2' and F('*CC(*)O') == 'C2H4O'
assert 0.5 * 72 == 36.0 and 0.5 * 44 == 22.0 and 5*12 + 8 + 32 == 100
assert not check_item('POL-RES-013')
print('POL-RES-013 OK')
