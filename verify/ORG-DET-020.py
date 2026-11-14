"""ORG-DET-020 検証：酸の強さによる分離と反応式"""
import sys, json
sys.path.insert(0, '/home/claude/chem/work/g09')
from g09lib import report_equations
pKa = {'benzoic': 4.20, 'carbonic': 6.35, 'phenol': 9.99}
def liberated_by(acid_pKa): return [k for k in ('benzoic', 'phenol') if pKa[k] > acid_pKa]
assert liberated_by(pKa['carbonic']) == ['phenol']           # CO2
assert liberated_by(-6) == ['benzoic', 'phenol']              # HCl
nahco3 = [k for k in ('benzoic', 'phenol') if pKa[k] < pKa['carbonic']]
assert nahco3 == ['benzoic']
T = {t['id']: set(t['pos']) for t in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/tests.json'))['tests']}
assert 'benzoic-acid' in T['nahco3'] and 'phenol' not in T['nahco3']
print('acid-strength logic OK')
report_equations('ORG-DET-020')
