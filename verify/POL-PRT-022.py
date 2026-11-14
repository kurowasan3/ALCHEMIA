# POL-PRT-022 検証：立体異性体の数
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from g12lib import check_item
def count(smi):
    m = Chem.MolFromSmiles(smi)
    s = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m, options=StereoEnumerationOptions(unique=True))}
    def mirror(x):
        return Chem.CanonSmiles(x.replace('@@', '!').replace('@', '@@').replace('!', '@'))
    return len(s), sum(mirror(x) == x for x in s)
r = {k: count(v) for k, v in {'鎖状Ala-Ala': 'NC(C)C(=O)NC(C)C(=O)O', '環状Ala-Ala': 'CC1NC(=O)C(C)NC1=O', '環状Ala-Ser': 'CC1NC(=O)C(CO)NC1=O'}.items()}
print(r)
assert r == {'鎖状Ala-Ala': (4, 0), '環状Ala-Ala': (3, 1), '環状Ala-Ser': (4, 0)}
assert not check_item('POL-PRT-022')
print('POL-PRT-022 OK')
