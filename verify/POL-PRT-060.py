# 検証スクリプト POL-PRT-060  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/POL-PRT-060/item.json'
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
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


from rdkit import Chem
from itertools import permutations, product
cs = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
FR = {'G': 'NCC(=O)', 'A': 'NC(C)C(=O)', 'Da': 'NC(CC(=O)O)C(=O)', 'Db': 'NC(C(=O)O)CC(=O)'}
def linear():
    out = {}
    for seq in permutations('GAD'):
        opts = [['Da','Db'] if (r=='D' and i<2) else (['Da'] if r=='D' else [r]) for i, r in enumerate(seq)]
        for combo in product(*opts):
            out.setdefault(cs(''.join(FR[c] for c in combo) + 'O'), combo)
    return out
L = linear(); print(len(L)); assert len(L) == 10
def cyclic():
    out = {}
    for seq in permutations('GAD'):
        for d in ('Da','Db'):
            combo = [d if r=='D' else r for r in seq]
            s = ''.join(FR[c] for c in combo)
            s = 'N9' + s[1:]
            assert s.endswith('C(=O)'); s = s[:-5] + 'C9=O'
            out.setdefault(cs(s), combo)
    return out
C = cyclic(); print(len(C)); assert len(C) == 4
def ring_size(s):
    m = Chem.MolFromSmiles(s); ri = m.GetRingInfo()
    return max(len(r) for r in ri.AtomRings())
sizes = sorted({ring_size(s) for s in C})
print(sizes); assert sizes == [9, 10]
X = [c for c in L.values() if c[0]=='G' and c[1] in ('Da','Db') and c[2]=='A']
print('X', X); assert len(X) == 2
Y = [c for c in L.values() if c[0]=='A' and c[2]=='Da']
print('Y', Y); assert Y == [('A','G','Da')]
# 示性式との照合
assert cs('NCC(=O)NC(CC(=O)O)C(=O)NC(C)C(=O)O') in L and cs('NCC(=O)NC(C(=O)O)CC(=O)NC(C)C(=O)O') in L
assert cs('NC(C)C(=O)NCC(=O)NC(C(=O)O)CC(=O)O') == cs(FR['A']+FR['G']+'NC(CC(=O)O)C(=O)O')
print('OK')
