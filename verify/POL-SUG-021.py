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
# ---- POL-SUG-021 アルドースの立体配置：RDKit でメソ判定し，割り当ての一意性を確認
from itertools import permutations
from rdkit import Chem
def chain(top, conf, bottom):
    # フィッシャー投影式を上から書く： [C@H](O) = OH 右（D-グリセルアルデヒド=R で較正）
    s = top
    for c in conf: s += '[C@H](O)' if c == 'R' else '[C@@H](O)'
    return s + bottom
assert Chem.FindMolChiralCenters(Chem.MolFromSmiles(chain('O=C', 'R', 'CO')))[0][1] == 'R'
# D-グルコース (R,L,R,R) の CIP = 2R,3S,4R,5R を確認
assert [c for _, c in Chem.FindMolChiralCenters(Chem.MolFromSmiles(chain('O=C', 'RLRR', 'CO')))] == ['R', 'S', 'R', 'R']
def achiral(s):
    a = Chem.MolToSmiles(Chem.MolFromSmiles(s))
    b = Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@', '!').replace('@', '@@').replace('!', '@')))
    return a == b
pent = ['RRR', 'LRR', 'RLR', 'LLR']   # C2,C3,C4
def exp1(p): return achiral(chain('OC', p, 'CO'))
def exp2(p): return sorted(achiral(chain('OC(=O)', c + p, 'C(=O)O')) for c in 'RL')
obs = {'P': (True, [False, True]), 'Q': (False, [False, False]), 'R': (True, [False, False]), 'S': (False, [False, True])}
sols = []
for perm in permutations(pent):
    asg = dict(zip('PQRS', perm))
    if all((exp1(asg[k]), exp2(asg[k])) == obs[k] for k in 'PQRS'): sols.append(asg)
print('解:', sols)
assert sols == [{'P': 'RRR', 'Q': 'LRR', 'R': 'RLR', 'S': 'LLR'}]
# (3) 片方の実験だけの場合
def amb(fn):
    groups = {}
    for p in pent: groups.setdefault(str(fn(p)), []).append(p)
    return sorted(v for v in groups.values())
print('実験1のみ', amb(exp1), '実験2のみ', amb(exp2))
assert amb(exp1) == [['LRR', 'LLR'], ['RRR', 'RLR']]
assert amb(exp2) == [['LRR', 'RLR'], ['RRR', 'LLR']]
# (2) P,S から光学不活性ジカルボン酸を与えるヘキソース
assert achiral(chain('OC(=O)', 'RRRR', 'C(=O)O')) and not achiral(chain('OC(=O)', 'LRRR', 'C(=O)O'))
assert achiral(chain('OC(=O)', 'RLLR', 'C(=O)O')) and not achiral(chain('OC(=O)', 'LLLR', 'C(=O)O'))
# (4) D-アルドヘキソース 8 種の還元体
from itertools import product
n = [''.join(t) + 'R' for t in product('RL', repeat=3) if achiral(chain('OC', ''.join(t) + 'R', 'CO'))]
print('メソのアルジトールを与える D-ヘキソース', n)
assert sorted(n) == ['RLLR', 'RRRR']
it = load('POL-SUG-021'); check_all_equations(it); check_subitems(it)
for s_ in ['P：2 位 右，3 位 右', 'Q：2 位 左，3 位 右', 'R：2 位 右，3 位 左', 'S：2 位 左，3 位 左', '(4) 2 種類']: in_text(s_, it)
print('POL-SUG-021 OK')
