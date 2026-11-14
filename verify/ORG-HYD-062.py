# 検証スクリプト ORG-HYD-062   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/ORG-HYD-062/item.json'
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
def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def check_structure():
    n_items = len(re.findall(r'\\item\b', re.sub(r'\\begin\{itemize\}.*?\\end\{itemize\}', '', it['problem'], flags=re.S)))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items + 1)), nums
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
import numpy as np, itertools
def nstereo(s): return len({Chem.MolToSmiles(m) for m in EnumerateStereoisomers(Chem.MolFromSmiles(s))})
assert nstereo('CC=CC=CC=CC') == 6      # (1) 2,4,6-オクタトリエンの立体異性体
assert nstereo('CC=CC=CC') == 3          # 2,4-ヘキサジエン
# 1,2-ジメチル環の 2 つのジアステレオマーについて，3D 構造からメチル基の向き（シス/トランス）と
# 鏡像体との重ね合わせ（メソ体かどうか）を判定する
def cis_and_achiral(smi, ring_atoms, me):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=7); AllChem.MMFFOptimizeMolecule(m)
    X = m.GetConformer().GetPositions()
    R = X[ring_atoms]; c = R.mean(0); n = np.linalg.svd(R - c)[2][2]
    s = [np.dot(X[i] - c, n) for i in me]
    cis = s[0] * s[1] > 0
    mol = Chem.MolFromSmiles(smi)
    mirror = Chem.MolToSmiles(Chem.MolFromSmiles(smi.replace('@@', '!').replace('@', '@@').replace('!', '@')))
    achiral = (Chem.MolToSmiles(mol) == mirror)
    return cis, achiral
for smi in ['C[C@H]1CCCC[C@@H]1C', 'C[C@H]1CCCC[C@H]1C']:
    cis, ach = cis_and_achiral(smi, [1, 2, 3, 4, 5, 6], [0, 7])
    print('1,2-ジメチルシクロヘキサン', smi, 'シス' if cis else 'トランス', 'メソ' if ach else 'キラル')
    assert cis == ach        # シス体 = メソ体，トランス体 = キラル
for smi in ['C[C@H]1C=CC=C[C@@H]1C', 'C[C@H]1C=CC=C[C@H]1C']:
    cis, ach = cis_and_achiral(smi, [1, 2, 3, 4, 5, 6], [0, 7])
    print('5,6-ジメチル-1,3-シクロヘキサジエン', smi, 'シス' if cis else 'トランス', 'メソ' if ach else 'キラル')
    assert cis == ach
# 水素化で 5,6 位の立体が保たれる（C=C だけに H2 を付加させても立体中心のタグは変わらない）
for smi in ['C[C@H]1C=CC=C[C@@H]1C', 'C[C@H]1C=CC=C[C@H]1C']:
    m = Chem.RWMol(Chem.MolFromSmiles(smi))
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE: b.SetBondType(Chem.BondType.SINGLE)
    Chem.SanitizeMol(m); h = Chem.MolToSmiles(m)
    print('水素化', smi, '→', h, '鏡像体をもたない' if Chem.MolToSmiles(Chem.MolFromSmiles(h.replace('@@','!').replace('@','@@').replace('!','@'))) == h else 'キラル')
# モデル（回転様式 × 外側/内側）による生成物の判定
def mode(npi, cond):
    return {(4, 'heat'): 'con', (4, 'light'): 'dis', (6, 'heat'): 'dis', (6, 'light'): 'con'}[(npi, cond)]
def ring_rel(pos1, pos2, md):      # pos: 'out' / 'in'
    same = (pos1 == pos2)
    return ('cis' if same else 'trans') if md == 'dis' else ('trans' if same else 'cis')
pos = {'t': 'out', 'c': 'in'}
res = {}
for a, b, c in itertools.product('tc', repeat=3):
    if b != 'c': continue
    key = tuple(sorted([a, c]))
    for cond in ('heat', 'light'):
        res[(key, cond)] = ring_rel(pos[a], pos[c], mode(6, cond))
print(res)
assert res[(('t', 't'), 'heat')] == 'cis' and res[(('t', 't'), 'light')] == 'trans'
assert res[(('c', 't'), 'heat')] == 'trans' and res[(('c', 't'), 'light')] == 'cis'
assert res[(('c', 'c'), 'heat')] == 'cis' and res[(('c', 'c'), 'light')] == 'trans'
# (3) 加熱→トランス（C はキラル），光→シス（E はメソ）を満たす A は一意
cand = [k for k in {('t', 't'), ('c', 't'), ('c', 'c')} if res[(k, 'heat')] == 'trans' and res[(k, 'light')] == 'cis']
print('A の候補', cand); assert cand == [('c', 't')]
# (4) シクロブテンの加熱開環（同旋的）
opens = {'cis': [], 'trans': []}
for p1, p2 in itertools.product(['out', 'in'], repeat=2):
    rel = ring_rel(p1, p2, mode(4, 'heat'))
    diene = tuple(sorted('t' if p == 'out' else 'c' for p in (p1, p2)))
    if diene not in opens[rel]: opens[rel].append(diene)
print('開環', opens)
assert opens['cis'] == [('c', 't')] and sorted(opens['trans']) == [('c', 'c'), ('t', 't')]
print('すべて OK')
