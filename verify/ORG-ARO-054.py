# 検証スクリプト ORG-ARO-054（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
    n_eq = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
    print('検査した式の数:', n_eq)
    for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem5/out/items/ORG-ARO-054/item.json')

from itertools import permutations, combinations
from rdkit import Chem
def hex_count(ms):
    G = [tuple((i + r) % 6 for i in range(6)) for r in range(6)] + [tuple((r - i) % 6 for i in range(6)) for r in range(6)]
    seen = set(); n = 0
    for p in set(permutations(ms)):
        if p in seen: continue
        n += 1
        for g in G: seen.add(tuple(p[g[i]] for i in range(6)))
    return n
res = [hex_count(list(x)) for x in ['XXXHHH', 'XXYHHH', 'XYZHHH', 'YXXXHH', 'XXXXHH', 'XXHHHH']]
print(res); assert res == [3, 6, 10, 6, 3, 3]
def rd_count(subs):
    out = set()
    for perm in set(permutations(subs + ['[H]'] * (6 - len(subs)))):
        smi = 'c1' + ''.join(f'c({s})' if i else '' for i, s in enumerate(perm))
        atoms = ['c(%s)' % s if s != '[H]' else 'c' for s in perm]
        smi = atoms[0][0] + '1' + atoms[0][1:] + ''.join(atoms[1:]) + '1'
        out.add(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
    return len(out)
r2 = [rd_count(['C', 'C', 'C']), rd_count(['C', 'C', 'O']), rd_count(['C', 'O', '[N+](=O)[O-]']), rd_count(['O'] + ['[N+](=O)[O-]'] * 3)]
print(r2); assert r2 == [3, 6, 10, 6]

print('ALL OK')
