import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

def load(item_id):
    return json.load(open(f'/home/claude/chem/out2/items/{item_id}/item.json', encoding='utf-8'))

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

def check_all_equations(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return len(eqs)

def check_subitems(it):
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    if n_items:
        assert [int(x) for x in nums] == list(range(1, n_items + 1)), (n_items, nums)
    print('小問数OK:', n_items)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def in_text(val_str, it, key='answer'):
    assert val_str in it[key], (val_str, key)
# ---- POL-SUG-020 三糖のメチル化分析：候補の網羅と絞り込み（抽象モデル）
from itertools import product, permutations
units = ['Gal', 'GlcA', 'GlcB']  # two glucoses are identical in kind
POS = [2, 3, 4, 5, 6]
def methyl_set(link_in, ring):
    return frozenset(p for p in POS if p not in link_in and p != ring)
def candidates():
    res = {}
    # trees over 3 labelled nodes with root (reducing end); parent maps
    for root in units:
        others = [u for u in units if u != root]
        shapes = []
        a, b = others
        shapes.append({a: root, b: root})       # branched at root
        shapes.append({a: root, b: a})          # chain root<-a<-b
        shapes.append({b: root, a: b})
        for par in shapes:
            for links in product(POS, repeat=2):
                lk = dict(zip(others, links))
                # parent positions must be distinct per parent
                inc = {u: [] for u in units}
                for ch, p in par.items(): inc[p].append(lk[ch])
                if any(len(v) != len(set(v)) for v in inc.values()): continue
                for rings in product([4, 5], repeat=2):
                    rg = dict(zip(others, rings))
                    if any(rg[u] in inc[u] for u in others): continue
                    # canonical form: replace GlcA/GlcB labels by Glc and sort
                    def node(u):
                        kids = sorted(node(c) + (lk[c],) for c in par if par[c] == u)
                        return (u[:3], rg.get(u, None), tuple(kids))
                    res.setdefault(node(root), []).append((node(root), root, dict(par), dict(lk), dict(rg)))
    return res
def exp2(c):
    tree, root, par, lk, rg = c
    inc = {u: [lk[ch] for ch, p in par.items() if p == u] for u in units}
    outs = set()
    for rr in [4, 5]:
        if rr in inc[root]: continue
        prods = []
        for u in units:
            ring = rg[u] if u != root else rr
            prods.append((u[:3], methyl_set(inc[u], ring)))
        outs.add(tuple(sorted(prods, key=lambda x: (x[0], sorted(x[1])))))
    return outs
def exp3(c):
    tree, root, par, lk, rg = c
    inc = {u: [lk[ch] for ch, p in par.items() if p == u] for u in units}
    prods = []
    for u in units:
        if u == root: prods.append((u[:3] + '-onic', methyl_set(inc[u], None)))
        else: prods.append((u[:3], methyl_set(inc[u], rg[u])))
    return tuple(sorted(prods, key=lambda x: (x[0], sorted(x[1]))))
def key(c): return c[0]
obs2 = tuple(sorted([('Gal', frozenset({2,3,4,6})), ('Glc', frozenset({2,3,4})), ('Glc', frozenset({2,3,6}))], key=lambda x: (x[0], sorted(x[1]))))
obs3 = tuple(sorted([('Gal', frozenset({2,3,4,6})), ('Glc', frozenset({2,3,4})), ('Glc-onic', frozenset({2,3,5,6}))], key=lambda x: (x[0], sorted(x[1]))))
cands = [v for vs in candidates().values() for v in vs]
ok2 = {}
for c in cands:
    if obs2 in exp2(c) and c[1] != 'Gal': ok2[key(c)] = c   # reducing end must be glucose? check below
ok2_all = {key(c): c for c in cands if obs2 in exp2(c)}
it = load('POL-SUG-020'); check_all_equations(it); check_subitems(it)
assert len(ok2_all) == 4, len(ok2_all)
names = {}
for k, c in ok2_all.items():
    print(k, '-> exp3', sorted((p, sorted(s)) for p, s in exp3(c)))
ok3 = [k for k, c in ok2_all.items() if exp3(c) == obs3]
assert len(ok3) == 1
# 期待：還元末端Glc <-4- Glc(六員環) <-6- Gal(六員環)
assert ok3[0] == ('Glc', None, (('Glc', 5, (('Gal', 5, (), 6),), 4),)), ok3
# 実験3で区別できない組 = (c),(d)
groups = {}
for k, c in ok2_all.items(): groups.setdefault(exp3(c), []).append(k)
dup = [v for v in groups.values() if len(v) > 1]
assert len(dup) == 1 and len(dup[0]) == 2
print('実験3で区別できない組:', dup[0])
# RDKit: T の分子式と OH の数
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
def unit(sugar='glc', anomer='a', subs=None, donor_tag=None):
    subs = subs or {}
    a1 = {'a': '@', 'b': '@@', None: ''}[anomer]
    a4 = {'glc': '@@', 'gal': '@'}[sugar]
    o = lambda p: 'O' + subs.get(p, '')
    first = donor_tag if donor_tag else o(1)
    c1 = f'[C{a1}H]' if a1 else '[CH]'
    return f"{first}{c1}1O[C@H](C{o(6)})[C{a4}H]({o(4)})[C@H]({o(3)})[C@H]1{o(2)}"
T = Chem.molzip(Chem.MolFromSmiles('.'.join([unit('gal', 'b', donor_tag='[*:1]'),
        unit('glc', 'a', subs={6: '[*:1]'}, donor_tag='[*:2]'), unit('glc', None, subs={4: '[*:2]'})])))
f = D.CalcMolFormula(T); mw = 18*12 + 32*1.0 + 16*16
nOH = len(T.GetSubstructMatches(Chem.MolFromSmarts('[OX2H]')))
print(Chem.MolToSmiles(T), f, mw, 'OH', nOH)
assert f == 'C18H32O16' and nOH == 11 and round(mw) == 504
print('POL-SUG-020 OK')
