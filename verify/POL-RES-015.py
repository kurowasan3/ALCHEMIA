# 検証スクリプト POL-RES-015  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/POL-RES-015/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, formula_mass
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it = json.load(open(ITEM, encoding='utf-8'))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 5, (n_items, n_ans)
eqs = sorted({b for k in ('problem', 'answer', 'solution') for b in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[k])
              if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm|\$', b)})
print(eqs); assert eqs == ['2C8H10O + CH2O -> C17H20O2 + H2O']
for e in eqs: assert check_equation(e) is None, check_equation(e)
def reactive_sites(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    c1 = o.GetNeighbors()[0].GetIdx()
    dm = Chem.GetDistanceMatrix(m)
    n = 0
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetSymbol() == 'C' and dm[c1][a.GetIdx()] in (1, 3) and any(nb.GetSymbol() == 'H' for nb in a.GetNeighbors()):
            n += 1
    return n
res = [reactive_sites(s) for s in ('Oc1ccccc1', 'Cc1ccc(O)cc1', 'Cc1ccc(O)c(C)c1')]; print('反応点', res)
assert res == [3, 2, 1]
for s, f in (('Oc1ccccc1', 'C6H6O'), ('Cc1ccc(O)cc1', 'C7H8O'), ('Cc1ccc(O)c(C)c1', 'C8H10O')):
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == f
dimer = 'Cc1cc(C)c(O)c(Cc2cc(C)cc(C)c2O)c1'
assert CalcMolFormula(Chem.MolFromSmiles(dimer)) == 'C17H20O2'
assert reactive_sites('Cc1cc(C)c(O)c(C)c1') == 0  # 結合した後の環にはもう反応点がない
M = formula_mass('C17H20O2', {'H': 1.0, 'C': 12, 'O': 16}); print('分子量', M); assert M == 256
print('全検証OK:', ITEM)
