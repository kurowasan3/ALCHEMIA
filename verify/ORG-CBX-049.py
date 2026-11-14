"""ORG-CBX-049 検証：2 候補の分子式，酸化的・還元的オゾン分解の生成物"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt

F1 = 'CCCC=CCCC=CCCCCCCC(=O)O'
F2 = 'CCCC=CCCCCCCC=CCCC(=O)O'
for f in (F1, F2):
    assert CalcMolFormula(Chem.MolFromSmiles(f)) == 'C16H28O2'

def cleave(smi, end):
    """全 C=C を切断し，切断端を end（'ox':COOH, 'red':CH2OH）にする"""
    m = Chem.RWMol(Chem.MolFromSmiles(smi))
    bonds = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE
             and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C']
    pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in bonds]
    for a, b in pairs:
        m.RemoveBond(a, b)
    for a, b in pairs:
        for x in (a, b):
            if end == 'ox':
                o1 = m.AddAtom(Chem.Atom(8)); o2 = m.AddAtom(Chem.Atom(8))
                m.AddBond(x, o1, Chem.BondType.DOUBLE); m.AddBond(x, o2, Chem.BondType.SINGLE)
            else:
                o = m.AddAtom(Chem.Atom(8)); m.AddBond(x, o, Chem.BondType.SINGLE)
    mm = m.GetMol()
    for a in mm.GetAtoms():
        a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    Chem.SanitizeMol(mm)
    return sorted(Chem.MolToSmiles(x) for x in Chem.GetMolFrags(mm, asMols=True))

def methylate(smis):
    rxn = AllChem.ReactionFromSmarts('[CX3:1](=[O:2])[OX2H1:3]>>[C:1](=[O:2])[O:3]C')
    out = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        while m.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OH]')):
            m = rxn.RunReactants((m,))[0][0]; Chem.SanitizeMol(m)
        out.append(Chem.MolToSmiles(m))
    return sorted(out)

ox1 = methylate(cleave(F1, 'ox')); ox2 = methylate(cleave(F2, 'ox'))
assert ox1 == ox2
exp1 = sorted(Chem.CanonSmiles(s) for s in ['CCCC(=O)OC', 'COC(=O)CCC(=O)OC', 'COC(=O)CCCCCCC(=O)OC'])
assert ox1 == exp1, ox1
red1 = cleave(methylate([F1])[0], 'red'); red2 = cleave(methylate([F2])[0], 'red')
diol4 = Chem.CanonSmiles('OCCCCO')
assert diol4 in red1 and diol4 not in red2
assert sorted(red1) == sorted(Chem.CanonSmiles(s) for s in ['CCCCO', 'OCCCCO', 'OCCCCCCCC(=O)OC'])
mw = 12 * 16 + 1.0 * 28 + 16 * 2
assert mw == 252 and abs(2.52 / mw * 2 * 22.4e3 - 448) < 1e-6
print('ORG-CBX-049 OK', red1, red2)
