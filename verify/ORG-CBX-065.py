# ORG-CBX-065 検証：AB2 型ヒドロキシ酸（2,2-ビス(ヒドロキシメチル)プロパン酸）のオリゴエステル
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, FindMolChiralCenters
from chemlib import check_equation
mass = {'C': 12, 'H': 1, 'O': 16}
def mw(m):
    m = Chem.AddHs(m); return sum(mass[a.GetSymbol()] for a in m.GetAtoms())
def trees(n):
    if n == 0: return [None]
    return [(a, b) for k in range(n) for a in trees(k) for b in trees(n - 1 - k)]
def core(t):
    arm = lambda c: 'CO' + ('C(=O)' + core(c) if c else '')
    return 'C(C)(' + arm(t[0]) + ')' + arm(t[1])
B = Chem.MolFromSmiles('CC(CO)(CO)C(=O)O')
assert rdMolDescriptors.CalcMolFormula(B) == 'C5H10O4' and mw(B) == 134
# 分子量 496 となる n
for n in range(1, 8):
    if 134 * n - 18 * (n - 1) + 32 - 18 == 496: N = n
assert N == 4
iso = {}
for t in trees(N):
    m = Chem.MolFromSmiles('COC(=O)' + core(t))
    iso[Chem.MolToSmiles(m)] = m
print('isomers', len(iso))
assert len(iso) == 3
patt = Chem.MolFromSmarts('[OX2H]')
for s, m in iso.items():
    oh = len(m.GetSubstructMatches(patt))
    ch = len(FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
    print(s, rdMolDescriptors.CalcMolFormula(m), mw(m), 'OH', oh, 'chiral', ch)
    assert mw(m) == 496 and oh == 5
chir = sorted(len(FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)) for m in iso.values())
assert chir == [1, 2, 3]
# どの n でも OH = n+1（n=1..6 の全異性体）
for n in range(1, 7):
    for t in trees(n):
        m = Chem.MolFromSmiles('COC(=O)' + core(t))
        assert len(m.GetSubstructMatches(patt)) == n + 1
# アセチル化：OH 1 個につき +42
mA = 4.96 / 496
print('acetylated mass', mA * (496 + 5 * 42))
assert abs(mA * (496 + 5 * 42) - 7.06) < 1e-9
# 不斉炭素 1 個の異性体は「メチルエステル単位 → 1 単位 → 末端 2 単位」
A = Chem.MolToSmiles(Chem.MolFromSmiles('COC(=O)C(C)(CO)COC(=O)C(C)(COC(=O)C(C)(CO)CO)COC(=O)C(C)(CO)CO'))
assert len(FindMolChiralCenters(iso[A], includeUnassigned=True, useLegacyImplementation=False)) == 1
# 加水分解の式（分子式で）
assert check_equation('C21H36O13 + 4H2O -> 4C5H10O4 + CH4O') is None
print('OK')
