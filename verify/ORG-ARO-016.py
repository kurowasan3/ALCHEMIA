"""ORG-ARO-016 検証：スチレン経路"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from rdkit import Chem
from rdkit.Chem import AllChem
from g09lib import report_equations, formula, iodoform, aldehyde, can, has
S = dict(A='CCc1ccccc1', B='C=Cc1ccccc1', D='CC(Cl)c1ccccc1', E='CC(O)c1ccccc1', F='OCCc1ccccc1', G='O=CCc1ccccc1')
for k, s in S.items(): print(k, formula(s))
assert formula(S['E']) == formula(S['F']) == 'C8H10O'
assert iodoform(S['E']) and not iodoform(S['F']) and has(S['F'], '[CH2][OH]') and aldehyde(S['G'])
DEH = AllChem.ReactionFromSmarts('[C;!H0:1][C:2][OH:3]>>[C:1]=[C:2].[O:3]')
for k in 'EF':
    prods = set()
    for ps in DEH.RunReactants((Chem.MolFromSmiles(S[k]),)):
        p = ps[0]; Chem.SanitizeMol(p); prods.add(Chem.MolToSmiles(p))
    print(k, 'dehydration ->', prods); assert prods == {can(S['B'])}
report_equations('ORG-ARO-016')
