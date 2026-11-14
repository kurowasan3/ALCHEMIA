"""ORG-CBX-014 検証：LiAlH4 還元・加水分解・18O の追跡"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from rdkit import Chem
from rdkit.Chem import AllChem
from g09lib import report_equations, formula, can, iodoform
RED_ESTER = AllChem.ReactionFromSmarts('[C:1](=[O:2])[O:3][C:4]>>[C:1][O:2].[O:3][C:4]')
RED_CO = AllChem.ReactionFromSmarts('[C:1]=[O:2]>>[C:1][O:2]')
HYD = AllChem.ReactionFromSmarts('[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]')
def run_all(smi, rxns):
    frags = [Chem.MolFromSmiles(smi)]
    changed = True
    while changed:
        changed = False
        for i, f in enumerate(frags):
            for rx in rxns:
                ps = rx.RunReactants((f,))
                if ps:
                    new = []
                    for p in ps[0]:
                        Chem.SanitizeMol(p); new.append(Chem.MolFromSmiles(Chem.MolToSmiles(p)))
                    frags = frags[:i] + new + frags[i + 1:]
                    changed = True; break
            if changed: break
    return sorted(Chem.MolToSmiles(f) for f in frags)
red = lambda s: run_all(s, [RED_ESTER, RED_CO])
hyd = lambda s: run_all(s, [HYD])
c = {'ア': 'CCC(=O)OCC', 'イ': 'CCCC=O', 'ウ': 'CC(=O)OCCC', 'エ': 'CCC(C)=O', 'オ': 'CCCC(=O)OC'}
R = {k: red(v) for k, v in c.items()}
for k in c: print(k, R[k])
same = [(a, b) for a in c for b in c if a < b and R[a] == R[b]]
print('same', same); assert same == [('ア', 'ウ')]
print('hydrolysis ア', hyd(c['ア']), 'ウ', hyd(c['ウ']))
assert can('CCC(=O)O') in hyd(c['ア']) and can('CC(=O)O') in hyd(c['ウ'])
X = {'diacetate': 'CC(=O)OCCOC(C)=O', 'oxalate': 'CCOC(=O)C(=O)OCC', 'glycolate': 'CC(=O)OCC(=O)OCC'}
for k, s in X.items():
    r, h = red(s), hyd(s)
    alcs = [x for x in set(h) if not Chem.MolFromSmiles(x).HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OH]'))]
    print(k, formula(s), 'red', sorted(set(r)), 'hyd', sorted(set(h)), 'alcohols', alcs, [iodoform(a) for a in alcs])
    assert formula(s) == 'C6H10O4' and set(r) == {can('CCO'), can('OCCO')}
ok = [k for k, s in X.items() if len(set(hyd(s))) == 2 and all(iodoform(a) for a in set(hyd(s))
      if not Chem.MolFromSmiles(a).HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OH]')))]
print('(4)', ok); assert ok == ['oxalate']
# 18O
lab = 'CCC(=O)[18O]CC'
r = red(lab); print('18O reduction', r)
assert any('18O' in x and formula(x.replace('[18OH]', 'O')) == 'C2H6O' for x in r)
report_equations('ORG-CBX-014')
