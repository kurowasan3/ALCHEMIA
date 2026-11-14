# 検証スクリプト THE-BND-050  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/THE-BND-050/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
check_structure(); check_all_equations()
# 供与 O（H 付き）と受容 O の最短距離（配座の中で最小）と環の員数
mols = {
 'o-ニトロフェノール': 'Oc1ccccc1[N+](=O)[O-]',
 'p-ニトロフェノール': 'Oc1ccc(cc1)[N+](=O)[O-]',
 'サリチル酸': 'Oc1ccccc1C(=O)O',
 'p-ヒドロキシ安息香酸': 'Oc1ccc(cc1)C(=O)O',
 'マレイン酸一価イオン': 'OC(=O)/C=C\\C(=O)[O-]',
 'フマル酸一価イオン': 'OC(=O)/C=C/C(=O)[O-]',
 'o-クレゾール': 'Oc1ccccc1C',
}
res = {}
for name, smi in mols.items():
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    cids = AllChem.EmbedMultipleConfs(m, 20, randomSeed=7); AllChem.MMFFOptimizeMoleculeConfs(m)
    donors = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='O' and any(n.GetSymbol()=='H' for n in a.GetNeighbors())]
    accs = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='O']
    best = (99, None)
    for d in donors:
        h = [n.GetIdx() for n in m.GetAtomWithIdx(d).GetNeighbors() if n.GetSymbol()=='H'][0]
        for ac in accs:
            if ac == d: continue
            ring = len(Chem.GetShortestPath(m, h, ac)) # H..acc 経路の原子数＝環の員数
            if ring < 5: continue  # 同じ COOH 内の O は除く
            for c in cids:
                p = m.GetConformer(c).GetPositions()
                dist = np.linalg.norm(p[d]-p[ac])
                if dist < best[0]: best = (dist, ring)
    res[name] = best; print(name, 'O…O 最短 %.2f Å' % best[0], '環', best[1])
assert res['o-ニトロフェノール'][1] == 6 and res['o-ニトロフェノール'][0] < 3.0
assert res['サリチル酸'][1] == 6 and res['サリチル酸'][0] < 3.0
assert res['マレイン酸一価イオン'][1] == 7 and res['マレイン酸一価イオン'][0] < 3.0
for k in ('p-ニトロフェノール','p-ヒドロキシ安息香酸'):
    assert res[k][0] > 5.0, k
# フマル酸の 2 つの COOH 間（異なる基の O どうし）
m = Chem.AddHs(Chem.MolFromSmiles('OC(=O)/C=C/C(=O)[O-]')); AllChem.EmbedMultipleConfs(m, 20, randomSeed=7); AllChem.MMFFOptimizeMoleculeConfs(m)
mm = Chem.MolFromSmiles('OC(=O)/C=C/C(=O)[O-]')
g1 = [0,2]; g2 = [6,7]
dmin = min(np.linalg.norm(c.GetPositions()[i]-c.GetPositions()[j]) for c in m.GetConformers() for i in g1 for j in g2)
print('フマル酸イオン 基間 O…O 最短 %.2f' % dmin); assert dmin > 3.3
# o-クレゾール：受容体は OH の O しかない
assert sum(1 for a in Chem.MolFromSmiles(mols['o-クレゾール']).GetAtoms() if a.GetSymbol()=='O') == 1
# pKa の比
pk = {'mal': (1.9, 6.2), 'fum': (3.0, 4.4)}
print('K1比', 10**(pk['fum'][0]-pk['mal'][0]), 'K2比', 10**(pk['mal'][1]-pk['fum'][1]))
print('THE-BND-050 OK')
