# ORG-STR-032 検証：糖アルコール（エリトリトール・キシリトール）の部分エステルの位置異性体の数（立体異性体は区別しない）
from itertools import product
from rdkit import Chem
def smi(n, groups):
    # n 個の炭素の直鎖ポリオール，各 OH を groups[i] で置換（'H','Ac','Pr'）
    sub = {'H': 'O', 'Ac': 'OC(=O)C', 'Pr': 'OC(=O)CC'}
    parts = ['C(' + sub[g] + ')' for g in groups]
    return Chem.CanonSmiles(''.join(parts), useChiral=0)
def count(n, filt):
    s = set()
    for g in product(['H', 'Ac', 'Pr'], repeat=n):
        if filt(g): s.add(smi(n, g))
    return len(s)
res = {}
for n, name in ((3, 'glycerol'), (4, 'erythritol'), (5, 'xylitol')):
    for k in range(1, n+1):
        res[(name, k)] = count(n, lambda g, k=k: g.count('Ac') == k and g.count('Pr') == 0)
res[('xylitol', 'Ac1Pr1')] = count(5, lambda g: g.count('Ac') == 1 and g.count('Pr') == 1)
res[('erythritol', 'Ac1Pr1')] = count(4, lambda g: g.count('Ac') == 1 and g.count('Pr') == 1)
print(res)
assert [res[('erythritol', k)] for k in range(1, 5)] == [2, 4, 2, 1]
assert [res[('xylitol', k)] for k in range(1, 6)] == [3, 6, 6, 3, 1]
assert res[('xylitol', 'Ac1Pr1')] == 10 and res[('erythritol', 'Ac1Pr1')] == 6
# 左右対称（自分自身と鏡映の位置関係が一致）なキシリトール三酢酸エステルの位置の組
sym = [c for c in __import__('itertools').combinations(range(1, 6), 3) if set(6 - i for i in c) == set(c)]
assert sym == [(1, 3, 5), (2, 3, 4)], sym
# 分子式：キシリトール C5H12O5 + 3 CH3COOH - 3 H2O = C11H18O8
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert CalcMolFormula(Chem.MolFromSmiles(smi(5, ['Ac', 'Ac', 'Ac', 'H', 'H']))) == 'C11H18O8'
print('ORG-STR-032 OK')
