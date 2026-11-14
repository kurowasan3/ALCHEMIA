# POL-RES-008 検証：高分子の分類
import sys, json, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
P = {p['id']: p for p in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/polymers.json'))['polymers']}
# 繰り返し単位（* は主鎖の続き）と 主鎖にヘテロ原子結合があるか
U = {'a': ('*CC(*)C', 'pp', False), 'b': ('*CC(*)Cl', 'pvc', False), 'c': ('*NCCCCCC(*)=O', 'nylon6', True),
     'd': ('*OCCOC(=O)c1ccc(cc1)C(*)=O', 'pet', True), 'e': (None, 'cellulose', True), 'f': (None, 'amylose', True),
     'g': ('*CC(*)OC(C)=O', 'pvac', False), 'h': ('*CC(*)C#N', 'acrylic', False), 'i': ('*CC(C)=CC*', 'natural-rubber', False),
     'j': (None, 'phenol-resin', False)}
form = {}
for k, (smi, pid, hetero) in U.items():
    f = P[pid].get('rf')
    if smi:
        f = re.sub(r'\*\d*', '', CalcMolFormula(Chem.MolFromSmiles(smi)))
        assert P[pid].get('rf') in (f, None), (k, f, P[pid]['rf'])
    form[k] = f
form['j'] = 'C7H6O'   # フェノール樹脂の構成元素 C,H,O（代表として）
els = {k: set(re.findall(r'[A-Z][a-z]?', f)) for k, f in form.items()}
print(form)
q1 = [k for k in U if els[k] == {'C', 'H'}]
q2 = [k for k in U if 'N' in els[k]]
q3 = [k for k in U if U[k][2]]
same = [(a, b) for a in U for b in U if a < b and form[a] == form[b]]
print(q1, q2, q3, same)
assert q1 == ['a', 'i'] and q2 == ['c', 'h'] and q3 == ['c', 'd', 'e', 'f'] and same == [('e', 'f')]
assert not check_item('POL-RES-008')
print('POL-RES-008 OK')
