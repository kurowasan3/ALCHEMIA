# 検証スクリプト ORG-STR-027  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-STR-027/item.json'
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
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from itertools import combinations
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
seen = {}
for s in ['CCCCC','CC(C)CC','CC(C)(C)C']:
    m = Chem.MolFromSmiles(s)
    for a, b in combinations(range(m.GetNumAtoms()), 2):
        rw = Chem.RWMol(m)
        for i in (a, b):
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(i, o, Chem.BondType.SINGLE)
        mm = rw.GetMol()
        try: Chem.SanitizeMol(mm)
        except Exception: continue
        seen[Chem.MolToSmiles(mm)] = mm
assert len(seen) == 14
def nchiral(m): return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
prim = [k for k, m in seen.items() if all(sum(1 for x in at.GetNeighbors()[0].GetNeighbors() if x.GetSymbol()=='C')==1 for at in m.GetAtoms() if at.GetSymbol()=='O')]
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
assert sorted(prim) == sorted(can(s) for s in ['OCCCCCO','OCC(C)CCO','CCC(CO)CO','CC(C)(CO)CO']), prim
both = []
for k in prim:
    acid = Chem.MolFromSmiles(Chem.MolToSmiles(seen[k]).replace('CO','C(=O)O') if False else k)
    rw = Chem.RWMol(Chem.MolFromSmiles(k))
    for at in list(rw.GetAtoms()):
        if at.GetSymbol()=='O':
            c = at.GetNeighbors()[0]
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c.GetIdx(), o, Chem.BondType.DOUBLE)
    am = rw.GetMol(); Chem.SanitizeMol(am)
    if nchiral(Chem.MolFromSmiles(k)) and nchiral(am): both.append((k, Chem.MolToSmiles(am)))
assert both == [(can('OCC(C)CCO'), can('OC(=O)C(C)CC(=O)O'))], both
two = {k: len(set(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m))) for k, m in seen.items() if nchiral(m) == 2}
assert two == {can('CC(O)C(O)CC'): 4, can('OCC(C)C(O)C'): 4, can('CC(O)CC(C)O'): 3}, two
print('異性体 14，第一級×2 が 4，酸化後も不斉 1，不斉 2 個の 3 種と立体異性体数 OK')

print('ALL OK')
