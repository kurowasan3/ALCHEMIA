# ORG-HYD-018 検証：ラジカル置換の選択性（H の数 × 1個あたりの反応性）による生成割合
import json, re, sys
from collections import defaultdict
from sympy import Rational as Q, N
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def mono(smi, rel, X='Cl'):
    """一置換体（構造異性体）ごとの重み {SMILES: 重み} と，各置換体の不斉炭素の有無"""
    m = Chem.MolFromSmiles(smi); w = defaultdict(lambda: Q(0))
    for a in m.GetAtoms():
        h = a.GetTotalNumHs()
        if not h: continue
        deg = a.GetDegree()   # 炭素の結合している炭素数＝級
        rw = Chem.RWMol(m); x = rw.AddAtom(Chem.Atom(17 if X == 'Cl' else 35)); rw.AddBond(a.GetIdx(), x, Chem.BondType.SINGLE)
        mm = rw.GetMol(); Chem.SanitizeMol(mm); w[Chem.MolToSmiles(mm)] += h * rel[deg]
    return dict(w)
def nst(s): return len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s))})
# (1) 反応性の比
r2 = (Q(72) / 4) / (Q(28) / 6); r3 = (Q(36) / 1) / (Q(64) / 9)
print('r2', N(r2, 4), 'r3', N(r3, 4)); assert round(float(r2), 1) == 3.9 and round(float(r3), 1) == 5.1
rel = {1: Q(1), 2: Q(39, 10), 3: Q(51, 10)}
# (2) C6H14 異性体ごとの一置換体の数と割合
hex = {'hexane': 'CCCCCC', '2-Me-pentane': 'CC(C)CCC', '3-Me-pentane': 'CCC(C)CC', '2,3-diMe-butane': 'CC(C)C(C)C', '2,2-diMe-butane': 'CC(C)(C)CC'}
for k, s in hex.items():
    w = mono(s, rel); tot = sum(w.values())
    print(k, len(w), {p: round(float(100 * v / tot), 1) for p, v in w.items()})
three = [k for k, s in hex.items() if len(mono(s, rel)) == 3]; assert sorted(three) == ['2,2-diMe-butane', 'hexane']
w = mono('CC(C)(C)CC', rel); tot = sum(w.values()); pc = sorted((round(float(100 * v / tot)) for v in w.values()), reverse=True)
assert pc == [45, 39, 15], pc
w2 = mono('CCCCCC', rel); tot2 = sum(w2.values()); assert sorted(round(float(100 * v / tot2)) for v in w2.values()) == [16, 42, 42]
# (3) X の一置換体のうち不斉炭素をもつもの → 立体異性体を含めて 4 種
assert sum(nst(p) for p in w) == 4
assert [p for p in w if nst(p) == 2] == [Chem.MolToSmiles(Chem.MolFromSmiles('CC(Cl)C(C)(C)C'))]
# (4) 2-メチルブタン
w = mono('CCC(C)C', rel); tot = sum(w.values()); pc = {p: round(float(100 * v / tot), 1) for p, v in w.items()}; print(pc)
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
assert pc[can('ClCC(C)CC')] == 27.4 and pc[can('CC(C)(Cl)CC')] == 23.3 and pc[can('CC(Cl)C(C)C')] == 35.6 and pc[can('ClCCC(C)C')] == 13.7
assert sum(nst(p) for p in w) == 6
# (5) 臭素化 1 : 82 : 1600
relB = {1: Q(1), 2: Q(82), 3: Q(1600)}
wb = mono('CCC(C)C', relB, 'Br'); tb = sum(wb.values()); t3 = wb[can('CC(C)(Br)CC')] / tb
print('Br tertiary', N(100 * t3, 4)); assert round(float(100 * t3)) == 90
it = json.load(open('/home/claude/chem/out2/items/ORG-HYD-018/item.json')); n = 0
for t in (it['problem'], it['answer'], it['solution']):
    for eq in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', t):
        if eq.count('->') == 1 and 'boxed' not in eq and 'R' not in eq:
            r = check_equation(eq); assert r is None, (eq, r); n += 1; print('eq ok', eq)
assert n >= 1
print('OK')
