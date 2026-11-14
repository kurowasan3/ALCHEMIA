# POL-PRT-059 検証：ジヌクレオシド一リン酸・短い鎖の異性体の数え上げ（総当たり＋正規化）
from itertools import product, permutations

def dimers(units, positions):
    """units: 2 個のヌクレオシドの名前．リン酸1個が両者の OH を1つずつ結ぶ．"""
    seen = set()
    for p, q in product(positions, positions):
        key = tuple(sorted([(units[0], p), (units[1], q)]))
        seen.add(key)
    return seen

ribo = ['2', '3', '5']
deoxy = ['3', '5']
AA = dimers(['A', 'A'], ribo); AU = dimers(['A', 'U'], ribo)
dAA = dimers(['A', 'A'], deoxy); dAU = dimers(['A', 'U'], deoxy)
print(len(AA), len(AU), len(dAA), len(dAU))
assert (len(AA), len(AU), len(dAA), len(dAU)) == (6, 9, 3, 4)
rna_type = [k for k in AU if sorted(p for _, p in k) == ['3', '5']]
print('RNA型', rna_type)
assert len(rna_type) == 2
# A2 分子どうし：RNA 型は 1 種類（ApA）
assert len([k for k in AA if sorted(p for _, p in k) == ['3', '5']]) == 1

def chains(multiset, link_types):
    """直鎖：上流の糖の (2 または 3) 位と下流の糖の 5 位をリン酸で結ぶ．向きがあるので逆順は別物．"""
    seen = set()
    for seq in set(permutations(multiset)):
        for links in product(link_types, repeat=len(seq) - 1):
            seen.add((seq, links))
    return seen

c1 = chains(['A', 'A', 'U'], ['3'])
c2 = chains(['A', 'A', 'U'], ['2', '3'])
c3 = chains(['A', 'G', 'U'], ['3'])
c4 = chains(['A', 'G', 'U'], ['2', '3'])
print(len(c1), len(c2), len(c3), len(c4))
assert (len(c1), len(c2), len(c3), len(c4)) == (3, 12, 6, 24)
# 公式との照合：同じヌクレオシド n 位置 → n(n+1)/2，異なる → n^2
for n in (2, 3):
    assert len(dimers(['A','A'], [str(i) for i in range(n)])) == n*(n+1)//2
    assert len(dimers(['A','U'], [str(i) for i in range(n)])) == n*n
# 分子式（ApU：アデノシン C10H13N5O4 + ウリジン C9H12N2O6 + H3PO4 - 2H2O）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
r = check_equation('C10H13N5O4 + C9H12N2O6 + H3PO4 -> C19H24N7O12P + 2H2O')
print(r); assert r is None
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
# ApU：A の 3' 位の O − P − U の 5' 位の O
ApU = 'Nc1ncnc2c1ncn2C1OC(CO)C(OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)C1O'
f = CalcMolFormula(Chem.MolFromSmiles(ApU)); print('ApU', f)
assert f == 'C19H24N7O12P'
print('OK')
