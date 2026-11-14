"""ORG-STR-042 検証：反応式のつり合い・不飽和度の変化・分子式の差"""
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def dou(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); c = {}
    for a in m.GetAtoms(): c[a.GetSymbol()] = c.get(a.GetSymbol(), 0) + 1
    return (2*c.get('C', 0) + 2 + c.get('N', 0) - c.get('H', 0) - c.get('Cl', 0)) / 2
R = {  # 反応: (有機反応物, 有機生成物)
 'a': (['C#C'], ['C=CCl']), 'b': (['C'], ['CCl']), 'c': (['CCO'], ['C=C']), 'd': (['CCO'], ['CCOCC']),
 'e': (['CCO'], ['CC=O']), 'f': (['CCOC(C)=O'], ['CC(=O)O', 'CCO']), 'g': (['c1ccccc1'], ['O=[N+]([O-])c1ccccc1']),
 'h': (['C=C'], ['CC'])}
dec = []
for k, (r, p) in R.items():
    dr = max(dou(s) for s in r); dp = max(dou(s) for s in p)
    print(k, [CalcMolFormula(Chem.MolFromSmiles(s)) for s in r], '->', [CalcMolFormula(Chem.MolFromSmiles(s)) for s in p], dr, '->', dp)
    if dp < dr: dec.append(k)
print('不飽和度が減少:', dec); assert dec == ['a', 'h']
# 分子式の差
def diff(a, b):
    A, _ = parse_formula(a); B, _ = parse_formula(b); return {e: A[e]-B.get(e, 0) for e in A if A[e]-B.get(e, 0)}
assert diff('C2H6O', 'C2H4') == {'H': 2, 'O': 1}; assert diff('C2H6O', 'C2H4O') == {'H': 2}
it = json.load(open('/home/claude/chem5/out/items/ORG-STR-042/item.json'))
n = 0
for key in ('problem', 'answer', 'solution'):
    for m in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[key]):
        if re.search(r'->|<=>', m) and '\\text' not in m:
            r = check_equation(m); print(m, r or 'OK'); assert r is None; n += 1
assert n == 8
print('ALL OK')
