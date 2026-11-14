"""ORG-ARO-014 検証：C8H10O 芳香族異性体の列挙と Cl 置換体の数"""
import sys
sys.path.insert(0, '/home/claude/chem/work/g09')
from g09lib import report_equations, formula, has, n_stereocenters, can
from arom import benzenes, ring_cl_isomers, n_ring_subs
subs = ['O', 'C', 'CC', 'OC', 'OCC', 'CO', 'CCO', 'C(C)O', 'COC']
allc = sorted(s for s in benzenes(subs) if formula(s) == 'C8H10O')
phen = [s for s in allc if has(s, 'c[OH]')]
alc = [s for s in allc if has(s, 'C[OH]')]
eth = [s for s in allc if not has(s, '[OH]')]
print(len(allc), 'phenols', len(phen), 'alcohols', len(alc), 'ethers', len(eth))
assert (len(allc), len(phen), len(alc), len(eth)) == (19, 9, 5, 5)
for s in phen: print(' ', s, 'subs', n_ring_subs(s), 'Cl isomers', len(ring_cl_isomers(s)))
two_di = [s for s in phen if n_ring_subs(s) == 2 and len(ring_cl_isomers(s)) == 2]
two_tri = [s for s in phen if n_ring_subs(s) == 3 and len(ring_cl_isomers(s)) == 2]
print('di', two_di, 'tri', two_tri)
assert two_di == [can('CCc1ccc(O)cc1')]
assert sorted(two_tri) == sorted([can('Cc1cccc(C)c1O'), can('Cc1cc(C)cc(O)c1')])
ch = [s for s in alc if n_stereocenters(s)]
print('chiral alcohols', ch); assert ch == [can('CC(O)c1ccccc1')]
report_equations('ORG-ARO-014')
