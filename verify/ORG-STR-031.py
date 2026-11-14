ITEM = '/home/claude/chem4/out/items/ORG-STR-031/item.json'
# 検証スクリプト（グループ g07）  python3 -I この.py
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
print('小問数', n_items, '答え', n_ans)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
# 問題の示性式 / 生徒の名称どおりの構造 / 正しい名称の構造
data = {
 'a': ('CC(CC)CC',      'CC(CC)CC',      'CCC(C)CC',     '3-メチルペンタン'),   # 2-エチルブタン
 'b': ('CCC(C)C',       'CC(C)CC',       'CC(C)CC',      '2-メチルブタン'),     # 3-メチルブタン = 2-メチルブタン
 'c': ('C=CCC',         'CCC=C',         'C=CCC',        '1-ブテン'),           # 3-ブテン
 'd': ('CC(C)C(C)C',    'CC(C)C(C)C',    'CC(C)C(C)C',   '2,3-ジメチルブタン'),
 'e': ('CC=CC(C)C',     'CC(C)C=CC',     'CC=CC(C)C',    '4-メチル-2-ペンテン'), # 2-メチル-3-ペンテン
 'f': ('CC(C)(C)CC',    'CCC(C)(C)C',    'CC(C)(C)CC',   '2,2-ジメチルブタン'),  # 3,3-ジメチルブタン
}
for k, (prob, student, correct, nm) in data.items():
    assert can(prob) == can(student) == can(correct), k
print('構造の一致 OK')

# 簡易命名：最長鎖（二重結合を含む）と位置番号の最小化
import itertools
NUM = {4: 'ブ', 5: 'ペン', 6: 'ヘキ'}
SUB = {1: 'メチル', 2: 'エチル'}
KAZU = {1: '', 2: 'ジ', 3: 'トリ'}
def name(smi):
    m = Chem.MolFromSmiles(smi)
    n = m.GetNumAtoms()
    dbl = [b for b in m.GetBonds() if b.GetBondTypeAsDouble() == 2]
    best = None
    for i, j in itertools.permutations(range(n), 2):
        for path in [list(Chem.GetShortestPath(m, i, j))]:
            if dbl and not all(any(path[k] == b.GetBeginAtomIdx() and path[k+1] == b.GetEndAtomIdx() or path[k] == b.GetEndAtomIdx() and path[k+1] == b.GetBeginAtomIdx() for k in range(len(path)-1)) for b in dbl):
                continue
            L = len(path)
            dpos = []
            for b in dbl:
                for k in range(L-1):
                    if {path[k], path[k+1]} == {b.GetBeginAtomIdx(), b.GetEndAtomIdx()}: dpos.append(k+1)
            subs = []
            for k, a in enumerate(path):
                for nb in m.GetAtomWithIdx(a).GetNeighbors():
                    if nb.GetIdx() not in path:
                        # 枝の炭素数
                        seen = {a}; stack = [nb.GetIdx()]; c = 0
                        while stack:
                            x = stack.pop()
                            if x in seen: continue
                            seen.add(x); c += 1
                            stack += [y.GetIdx() for y in m.GetAtomWithIdx(x).GetNeighbors()]
                        subs.append((k+1, c))
            key = (-L, sorted(dpos), sorted(p for p, _ in subs))
            if best is None or key < best[0]:
                best = (key, L, dpos, subs)
    _, L, dpos, subs = best
    pre = ''
    groups = {}
    for p, c in subs: groups.setdefault(c, []).append(p)
    parts = []
    for c in sorted(groups, key=lambda c: SUB[c]):
        ps = sorted(groups[c]); parts.append(','.join(map(str, ps)) + '-' + KAZU[len(ps)] + SUB[c])
    pre = '-'.join(parts)
    if dpos:
        base = NUM[L] + 'テン' if L != 4 else 'ブテン'
        base = {4: 'ブテン', 5: 'ペンテン', 6: 'ヘキセン'}[L]
        core = f'{dpos[0]}-{base}'
        nm = (pre + '-' if pre else '') + core
        if not pre and L == 4 and False: pass
    else:
        nm = pre + {4: 'ブタン', 5: 'ペンタン', 6: 'ヘキサン'}[L]
    return nm
for k, (prob, student, correct, nm) in data.items():
    got = name(prob)
    print(k, got, nm)
    assert got == nm, (k, got, nm)
g3 = name('CCC(C)=C'); print('(3)', g3); assert g3 == '2-メチル-1-ブテン'
F = {k: CalcMolFormula(Chem.MolFromSmiles(v[0])) for k, v in data.items()}
print(F)
assert sorted(k for k in F if F[k] == 'C6H14') == ['a', 'd', 'f']
print('すべて OK')
