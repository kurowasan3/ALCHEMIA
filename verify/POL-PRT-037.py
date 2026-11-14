# 検証スクリプト POL-PRT-037  python3 -I この.py
ITEM = '/home/claude/chem4/outh/items/POL-PRT-037/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
AAS = {'Cys':'NC(CS)C(=O)O','Lys':'NCCCCC(N)C(=O)O','Ala':'CC(N)C(=O)O','Gly':'NCC(=O)O','Ser':'NC(CO)C(=O)O'}
AA = {'Cys':121,'Lys':146,'Ala':89,'Gly':75,'Ser':105}
for k,s in AAS.items():
    assert round(MolWt(Chem.MolFromSmiles(s))) == AA[k], k
pep = lambda seq: sum(AA[a] for a in seq) - 18*(len(seq)-1)
X = 'Cys Lys Ala Cys Lys Gly Gly Cys Lys Ser Cys'.split()
assert pep(X) == 1086 and pep(X) - 1082 == 4
def matchings(items):
    if not items: return [[]]
    a = items[0]; out = []
    for i in range(1, len(items)):
        rest = items[1:i] + items[i+1:]
        for m in matchings(rest): out.append([(a, items[i])] + m)
    return out
assert len(matchings([1,4,8,11])) == 3 and len(matchings(list(range(6)))) == 15
def frags(seq):
    cuts = [i+1 for i,a in enumerate(seq[:-1]) if a == 'Lys']
    b = [0] + cuts + [len(seq)]
    return [list(range(b[k], b[k+1])) for k in range(len(b)-1)]
def products(seq, pairing):
    fs = frags(seq); parent = list(range(len(fs)))
    def find(x):
        while parent[x] != x: x = parent[x]
        return x
    where = {p: k for k, f in enumerate(fs) for p in f}
    nss = {}
    for a, b in pairing:
        fa, fb = find(where[a-1]), find(where[b-1]); 
        if fa != fb: parent[fa] = fb
    groups = {}
    for k in range(len(fs)): groups.setdefault(find(k), []).append(k)
    res = []
    for g in groups.values():
        m = sum(pep([seq[p] for p in fs[k]]) for k in g)
        ss = sum(1 for a, b in pairing if where[a-1] in g)
        res.append(m - 2*ss)
    return sorted(res)
assert [pep([X[p] for p in f]) for f in frags(X)] == [249, 320, 363, 208]
res = {tuple(map(tuple, m)): products(X, m) for m in matchings([1,4,8,11])}
print(res)
hit = [m for m, r in res.items() if r == [526, 610]]
assert hit == [((1, 8), (4, 11))]
assert len(set(map(tuple, res.values()))) == 3
Y = X.copy(); Y[4] = 'Gly'
resY = {tuple(map(tuple, m)): products(Y, m) for m in matchings([1,4,8,11])}
print(resY)
assert resY[((1,4),(8,11))] == resY[((1,8),(4,11))] == [1047]
assert resY[((1,11),(4,8))] == [455, 592]
print('ジスルフィド結合の組み合わせ OK')

print('ALL OK')
