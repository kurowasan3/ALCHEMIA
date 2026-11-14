"""ORG-ARO-015 検証：フェノール性 OH の有無と反応式"""
import sys, json
sys.path.insert(0, '/home/claude/chem/work/g09')
from g09lib import report_equations, has
c = {'ア': 'OC(=O)c1ccccc1O', 'イ': 'COC(=O)c1ccccc1O', 'ウ': 'CC(=O)Oc1ccccc1C(=O)O',
     'エ': 'O=C(Oc1ccccc1)c1ccccc1', 'オ': 'COC(=O)c1ccc(O)cc1', 'カ': 'COc1ccccc1'}
pos = [k for k, s in c.items() if has(s, 'c[OX2H1]')]
print('phenolic OH', pos); assert pos == ['ア', 'イ', 'オ']
T = {t['id']: set(t['pos']) for t in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/tests.json'))['tests']}
assert {'salicylic-acid', 'methyl-salicylate'} <= T['fecl3'] and 'aspirin' not in T['fecl3']
# エステル 1 mol あたりの NaOH：酸 1 + 生じたフェノール 1
n_phenyl = 1 + int(has(c['エ'], 'C(=O)Oc')); print('NaOH for phenyl benzoate', n_phenyl); assert n_phenyl == 2
report_equations('ORG-ARO-015')
