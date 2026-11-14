# 検証スクリプト ORG-ARO-050  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-ARO-050/item.json'
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
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ORG-ARO-050：環化三量化のモデルを全列挙 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import product, permutations
def benz(subs):
    parts = [('c1' if i == 0 else 'c') + (f'({g})' if g else '') for i, g in enumerate(subs)]
    smi = ''.join(parts); k = smi.rfind('c'); smi = smi[:k+1] + '1' + smi[k+1:]
    return Chem.CanonSmiles(smi)
def trimer(alk):
    out = set()
    for perm in set(permutations(alk)):
        for fl in product([0, 1], repeat=3):
            subs = []
            for (r1, r2), f in zip(perm, fl): subs += [r1, r2] if f == 0 else [r2, r1]
            out.add(benz(subs))
    return out
def mix(*alks):   # 混合物：3 分子の選び方すべて
    from itertools import combinations_with_replacement
    s = set()
    for c in combinations_with_replacement(alks, 3): s |= trimer(list(c))
    return s
can = Chem.CanonSmiles
P, A, B = ('C', ''), ('', ''), ('C', 'C')
# (1)
r1 = trimer([P, P, P]); print(r1)
assert r1 == {can('Cc1cc(C)cc(C)c1'), can('Cc1ccc(C)c(C)c1')}
assert can('Cc1cccc(C)c1C') not in r1
# (2)
r2 = mix(A, B); print(r2)
assert r2 == {can('c1ccccc1'), can('Cc1ccccc1C'), can('Cc1ccc(C)c(C)c1C'), can('Cc1c(C)c(C)c(C)c(C)c1C')}
# (3) アセチレン＋プロピンで生じるキシレン
r3 = mix(A, P); xyl = {s for s in r3 if CalcMolFormula(Chem.MolFromSmiles(s)) == 'C8H10'}; print(xyl)
assert xyl == {can('Cc1ccccc1C'), can('Cc1cccc(C)c1'), can('Cc1ccc(C)cc1')}
# (4) C6H10 のアルキン 7 種（RDKit で列挙して確認）
alk6 = {'C#CCCCC':('CCCC',''), 'CC#CCCC':('C','CCC'), 'CCC#CCC':('CC','CC'), 'C#CC(C)CC':('C(C)CC',''),
        'C#CCC(C)C':('CC(C)C',''), 'C#CC(C)(C)C':('C(C)(C)C',''), 'CC#CC(C)C':('C','C(C)C')}
for s in alk6: assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C6H10'
# C6H10 の鎖式アルキンが 7 種で全部であることを，炭素骨格に三重結合を置く全列挙で確認
from rdkit.Chem import AllChem
skel = ['CCCCCC', 'CC(C)CCC', 'CCC(C)CC', 'CC(C)(C)CC', 'CC(C)C(C)C']
found = set()
for sk in skel:
    m = Chem.MolFromSmiles(sk)
    for b in m.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetDegree() <= 2 and a2.GetDegree() <= 2:
            rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.TRIPLE)
            found.add(Chem.MolToSmiles(rw))
print('C6H10 アルキン', len(found)); assert len(found) == 7 and found == {can(s) for s in alk6}
counts = {s: len(trimer([v]*3)) for s, v in alk6.items()}; print(counts)
assert [s for s, n in counts.items() if n == 1] == ['CCC#CCC'] and all(n == 2 for s, n in counts.items() if s != 'CCC#CCC')
# (5) C5H8 のアルキン
alk5 = {'C#CCCC':('CCC',''), 'CC#CCC':('C','CC'), 'C#CC(C)C':('C(C)C','')}
for s, v in alk5.items():
    r = trimer([v]*3)
    ringH = [Chem.MolFromSmiles(x) for x in r]
    noH = all(all(a.GetTotalNumHs() == 0 for a in m.GetAtoms() if a.GetIsAromatic()) for m in ringH)
    print(s, len(r), noH)
    assert (len(r) == 2 and noH) == (s == 'CC#CCC')
r5 = trimer([('C','CC')]*3)
assert r5 == {can('CCc1c(C)c(CC)c(C)c(CC)c1C'), can('CCc1c(C)c(C)c(CC)c(CC)c1C')}
for x in r5: assert CalcMolFormula(Chem.MolFromSmiles(x)) == 'C15H24'
print('ALL OK')
