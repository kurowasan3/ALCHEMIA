# POL-PRT-043 検証：Glu・Cys・Gly 各1個の鎖状トリペプチドの候補の網羅と，実験による絞り込み
import itertools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# 残基の官能基: N（アミノ基）と C（カルボキシ基，α か側鎖 γ か）
RES = {
 'Glu': ('NC(CCC(=O)O)C(=O)O', {'N': {'a': 0}, 'C': {'g': (4, 6), 'a': (7, 9)}}),
 'Cys': ('NC(CS)C(=O)O', {'N': {'a': 0}, 'C': {'a': (4, 6)}}),
 'Gly': ('NCC(=O)O', {'N': {'a': 0}, 'C': {'a': (2, 4)}}),
}
names = ['Glu', 'Cys', 'Gly']
combo = Chem.MolFromSmiles('.'.join(RES[n][0] for n in names))
off = [0, 10, 17]
cg = [(i, k, off[i] + c, off[i] + oh) for i, n in enumerate(names) for k, (c, oh) in RES[n][1]['C'].items()]
ng = [(i, k, off[i] + a) for i, n in enumerate(names) for k, a in RES[n][1]['N'].items()]
edges = [(c, n) for c in cg for n in ng if c[0] != n[0]]
cands = {}
for es in itertools.combinations(edges, 2):
    if len({e[0][:2] for e in es}) < 2 or len({e[1][:2] for e in es}) < 2: continue
    if {es[0][0][0], es[0][1][0]} == {es[1][0][0], es[1][1][0]}: continue   # 同じ2残基を2本で結ぶ（環）は除く → 3残基の木
    rw = Chem.RWMol(combo)
    for (ci, ck, c, oh), (ni, nk, n) in es: rw.AddBond(c, n, Chem.BondType.SINGLE)
    for a in sorted([e[0][3] for e in es], reverse=True): rw.RemoveAtom(a)
    m = rw.GetMol(); Chem.SanitizeMol(m)
    s = Chem.MolToSmiles(m)
    bonds = tuple(sorted((names[e[0][0]] + '-' + e[0][1], names[e[1][0]]) for e in es))
    cands[s] = bonds
print('全候補', len(cands))
for s, b in cands.items(): print(' ', b, s)
assert len(cands) == 12
for s in cands: assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C10H17N3O6S'

# (b) N 末端標識で Glu が標識される = Glu の N が結合に使われていない
b = {s: bd for s, bd in cands.items() if not any(acc == 'Glu' for _, acc in bd)}
print('(b) 後', len(b)); assert len(b) == 6
# (c) 酵素P: α-カルボキシ基によるペプチド結合だけを切る. 生成物が「遊離の Gly」と P で切れないジペプチド
def products_after_P(bd):
    keep = [x for x in bd if x[0] == 'Glu-g']   # 切れずに残る結合
    # 残る結合で連結した残基の組
    groups = [{n} for n in names]
    for donor, acc in keep:
        d = donor.split('-')[0]
        g1 = next(g for g in groups if d in g); g2 = next(g for g in groups if acc in g)
        if g1 is not g2: g1 |= g2; groups.remove(g2)
    return sorted(tuple(sorted(g)) for g in groups)
c = {s: bd for s, bd in b.items() if products_after_P(bd) == sorted([('Gly',), ('Cys', 'Glu')])}
print('(c) 後', len(c))
for s, bd in c.items(): print('  ', bd, s)
assert len(c) == 2
# (c) で除かれた候補の生成物も表示
for s, bd in b.items(): print('   P生成物', bd, products_after_P(bd))
# (d) 試薬Z: 同じ炭素に遊離の NH2 と遊離の COOH がある部分（α-アミノ酸部分）
patt = Chem.MolFromSmarts('[NX3;H2][CX4][CX3](=O)[OX2H1]')
d = {s: bd for s, bd in c.items() if Chem.MolFromSmiles(s).HasSubstructMatch(patt)}
print('(d) 後', d); assert len(d) == 1
ans = list(d)[0]
assert set(d[ans]) == {('Glu-g', 'Cys'), ('Cys-a', 'Gly')}   # γ-Glu-Cys-Gly（グルタチオン）
# 遊離の COOH は全候補で 2 個（中和滴定は区別に使えない）
for s in cands:
    assert len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts('C(=O)[OX2H1]'))) == 2
# 不斉炭素 2 個 → 立体異性体 4
m = Chem.MolFromSmiles(ans)
cc = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
print('不斉炭素', cc); assert len(cc) == 2
print('答えの構造', ans, rdMolDescriptors.CalcMolFormula(m))
print('ALL OK')
