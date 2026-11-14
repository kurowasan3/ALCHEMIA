# 検証スクリプト ORG-STR-053  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/ORG-STR-053/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0:
            return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                n += 1
    print('反応式', n, '本 OK')

def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items:
        assert n_items == n_ans, (n_items, n_ans)
    else:
        assert n_ans == 0
    print('小問数', n_items, '答え数', n_ans)

def check_svgs():
    import os
    d = os.path.dirname(ITEM)
    for k in ('problem_fig', 'solution_fig'):
        f = it[k]
        if f:
            s = open(os.path.join(d, f), encoding='utf-8').read()
            assert 'viewBox' in s and '<script' not in s and 'foreignObject' not in s
            assert not re.search(r'\son[a-z]+=', s) and 'href="http' not in s and 'base64' not in s
    print('SVG OK')

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

check_structure(); check_equations(); check_svgs()

# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T
E = {frozenset(['Me']): 3.8, frozenset(['Me', 'OH']): 1.5, frozenset(['OH']): 1.0}
front = {'OH': 0, 'Me': 120, 'H': 240}
backs = {'meso': {'OH': 180, 'Me': 300, 'H': 60}, 'chiral': {'OH': 180, 'Me': 60, 'H': 300}}
res = {}
for name, b in backs.items():
    out = []
    for rot in (0, 120, 240):   # a, b, c の順（OH を q→r→p）
        bb = {k: (v + rot) % 360 for k, v in b.items()}
        e = w = 0.0
        for f, fa in front.items():
            for k, ka in bb.items():
                if 'H' in (f, k): continue
                d = min((fa - ka) % 360, (ka - fa) % 360)
                if d == 60:
                    e += E[frozenset([f, k])]; w += E[frozenset([f, k])]
                    if f == k == 'OH': e += -10.0
        meme = min((120 - bb['Me']) % 360, (bb['Me'] - 120) % 360)
        out.append((round(e, 1), round(w, 1), bb, meme))
    res[name] = out
    print(name, out)
assert [x[0] for x in res['meso']] == [3.0, -3.7, -3.7]
assert [x[1] for x in res['meso']] == [3.0, 6.3, 6.3]
assert [x[0] for x in res['chiral']] == [6.8, -5.2, -6.0]
assert [x[1] for x in res['chiral']] == [6.8, 4.8, 4.0]
# 最安定配座の CH3 の関係
m = min(res['meso'], key=lambda x: x[0]); assert m[3] == 60
c = min(res['chiral'], key=lambda x: x[0]); assert c[3] == 180
# (3) 差
assert round(3.0 - (-3.7), 1) == 6.7 and round(6.8 - (-6.0), 1) == 12.8
# 水中の最安定の OH の関係
wm = min(res['meso'], key=lambda x: x[1]); assert min((0 - wm[2]['OH']) % 360, wm[2]['OH'] % 360) == 180
wc = min(res['chiral'], key=lambda x: x[1]); assert min((0 - wc[2]['OH']) % 360, wc[2]['OH'] % 360) == 60
# RDKit: メソ体（R,S）は OCCO=180 で CCCC≈180，キラル体（R,R）は ≈±60
for smi, want in [('C[C@@H](O)[C@@H](O)C', 'meso'), ('C[C@@H](O)[C@H](O)C', 'chiral')]:
    mol = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(mol, randomSeed=7)
    cip = sorted(x[1] for x in Chem.FindMolChiralCenters(mol))
    assert (cip == ['R', 'S']) == (want == 'meso')
    conf = mol.GetConformer(); T.SetDihedralDeg(conf, 2, 1, 3, 4, 180.0)
    cccc = abs(T.GetDihedralDeg(conf, 0, 1, 3, 5))
    print(want, cip, round(cccc))
    assert (cccc > 150) == (want == 'meso')
# アセタール：シス体はメソ，トランス体はキラル
def mirror(m):
    s = Chem.MolToSmiles(m); return Chem.CanonSmiles(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
def same_as_mirror(m): return Chem.MolToSmiles(m) == mirror(m)
rw = Chem.RWMol(Chem.MolFromSmiles('C[C@@H](O)[C@@H](O)C'))
c0 = rw.AddAtom(Chem.Atom(6)); m1 = rw.AddAtom(Chem.Atom(6)); m2 = rw.AddAtom(Chem.Atom(6))
rw.AddBond(c0, m1, Chem.BondType.SINGLE); rw.AddBond(c0, m2, Chem.BondType.SINGLE)
rw.AddBond(c0, 2, Chem.BondType.SINGLE); rw.AddBond(c0, 4, Chem.BondType.SINGLE)
ac = rw.GetMol(); Chem.SanitizeMol(ac); Chem.AssignStereochemistry(ac, cleanIt=True, force=True)
print('メソ体からのアセタール', Chem.MolToSmiles(ac), same_as_mirror(ac))
assert same_as_mirror(ac)
rw = Chem.RWMol(Chem.MolFromSmiles('C[C@@H](O)[C@H](O)C'))
c0 = rw.AddAtom(Chem.Atom(6)); m1 = rw.AddAtom(Chem.Atom(6)); m2 = rw.AddAtom(Chem.Atom(6))
rw.AddBond(c0, m1, Chem.BondType.SINGLE); rw.AddBond(c0, m2, Chem.BondType.SINGLE)
rw.AddBond(c0, 2, Chem.BondType.SINGLE); rw.AddBond(c0, 4, Chem.BondType.SINGLE)
ac2 = rw.GetMol(); Chem.SanitizeMol(ac2)
print('キラル体からのアセタール', Chem.MolToSmiles(ac2), same_as_mirror(ac2))
assert not same_as_mirror(ac2)

print('ALL OK')
