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
# ---- POL-SUG-023 二糖の異性体の数え上げ（RDKit molzip で構造を作り正規化 SMILES で重複除去）
from rdkit import Chem
def unit(sugar='glc', anomer='a', subs=None, donor_tag=None):
    subs = subs or {}
    a1 = {'a': '@', 'b': '@@', None: ''}[anomer]
    a4 = {'glc': '@@', 'gal': '@'}[sugar]
    o = lambda p: 'O' + subs.get(p, '')
    first = donor_tag if donor_tag else (subs[1] + 'O' if 1 in subs else 'O')
    c1 = f'[C{a1}H]' if a1 else '[CH]'
    return f"{first}{c1}1O[C@H](C{o(6)})[C{a4}H]({o(4)})[C@H]({o(3)})[C@H]1{o(2)}"
def disacc(s1, s2):
    """s1 が供与体（1 位で結合），s2 が受容体．全組合せの正規化 SMILES の集合"""
    out = {}
    for an1 in 'ab':
        for pos in [1, 2, 3, 4, 6]:
            an2s = 'ab' if pos == 1 else [None]
            for an2 in an2s:
                smi = '.'.join([unit(s1, an1, donor_tag='[*:1]'), unit(s2, an2, subs={pos: '[*:1]'})])
                m = Chem.molzip(Chem.MolFromSmiles(smi))
                can = Chem.MolToSmiles(m)
                out.setdefault(can, []).append((s1, an1, pos, s2, an2))
    return out
def reducing(m):
    # 還元末端：C1 に OH とエーテル O(環) が結合した炭素（ヘミアセタール）があるか
    patt = Chem.MolFromSmarts('[CX4;R]([OX2H])[OX2;R]')
    return m.HasSubstructMatch(patt)
def summary(pairs):
    allc = {}
    for a, b in pairs:
        for k, v in disacc(a, b).items(): allc.setdefault(k, []).extend(v)
    nred = sum(1 for k in allc if reducing(Chem.MolFromSmiles(k)))
    return allc, nred
it = load('POL-SUG-023'); check_all_equations(it); check_subitems(it)
gg, r = summary([('glc', 'glc')])
print('Glc-Glc', len(gg), 'reducing', r)
assert len(gg) == 11 and r == 8
same = [v for v in gg.values() if len(v) > 1]
assert same == [[('glc', 'a', 1, 'glc', 'b'), ('glc', 'b', 1, 'glc', 'a')]] or len(same) == 1
# (2) 単位が区別できない＝ C2 対称：αα, ββ（原子の対称クラスで確認）
def symmetric(can):
    m = Chem.MolFromSmiles(can)
    ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    anomeric = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'C' and sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'O') == 2]
    return len(anomeric) == 2 and ranks[anomeric[0]] == ranks[anomeric[1]]
sym = [v[0] for k, v in gg.items() if v[0][2] == 1 and symmetric(k)]
print('対称な 1,1-二糖', sym)
assert sorted((x[1], x[4]) for x in sym) == [('a', 'a'), ('b', 'b')]
# (3) メチル化分析のグループ：生成物 = 結合位置で決まる
groups = {}
for k, v in gg.items():
    s1, an1, pos, s2, an2 = v[0]
    groups.setdefault(pos, []).append(k)
sizes = sorted((p, len(v)) for p, v in groups.items())
print('グループ', sizes)
assert sizes == [(1, 3), (2, 2), (3, 2), (4, 2), (6, 2)]
gx, r2 = summary([('glc', 'gal'), ('gal', 'glc')])
print('Glc/Gal', len(gx), 'reducing', r2)
assert len(gx) == 20 and r2 == 16
for s_ in ['(1) 11 種類，還元性を示すもの 8 種類', '(3) 5 グループ', '(4) 20 種類，還元性を示すもの 16 種類']: in_text(s_, it)
print('POL-SUG-023 OK')
