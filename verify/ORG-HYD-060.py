# 検証スクリプト ORG-HYD-060   python3 -I この.py
ITEM = '/home/claude/chemk/out/items/ORG-HYD-060/item.json'
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
def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def check_structure():
    n_items = len(re.findall(r'\\item\b', re.sub(r'\\begin\{itemize\}.*?\\end\{itemize\}', '', it['problem'], flags=re.S)))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items + 1)), nums
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
check_all_equations()
# ---- 個別の検証（RDKit）----
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem import rdMolDescriptors as D
def canon(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s))
def nstereo(s): return len({Chem.MolToSmiles(m) for m in EnumerateStereoisomers(Chem.MolFromSmiles(s))})
def nchiral(s):
    return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
def has_ez(s):
    m = Chem.MolFromSmiles(s)
    from rdkit.Chem import FindPotentialStereo
    return any(si.type == Chem.StereoType.Bond_Double for si in Chem.FindPotentialStereo(m))
def adducts(diene_smi):
    """反応 SMARTS で 1,2-付加と（共役なら）1,4-付加を機械的に生成"""
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(diene_smi)
    prods = set()
    r12 = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1](Br)[C:2]Br')
    r14 = AllChem.ReactionFromSmarts('[C:1]=[C:2]-[C:3]=[C:4]>>[C:1](Br)[C:2]=[C:3][C:4]Br')
    for r in (r12, r14):
        for ps in r.RunReactants((m,)):
            p = ps[0]
            try:
                Chem.SanitizeMol(p); prods.add(Chem.MolToSmiles(p))
            except Exception: pass
    return prods
iso = adducts('C=C(C)C=C')
print('イソプレンの付加体', iso)
expect = {canon('BrCC(C)(Br)C=C'), canon('C=C(C)C(Br)CBr'), canon('BrCC(C)=CCBr')}
assert iso == expect
for s in iso: assert D.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C5H8Br2'
# (2) 不斉炭素・シストランス・立体異性体の合計
assert nchiral('BrCC(C)(Br)C=C') == 1 and not has_ez('BrCC(C)(Br)C=C')
assert nchiral('C=C(C)C(Br)CBr') == 1 and not has_ez('C=C(C)C(Br)CBr')
assert nchiral('BrCC(C)=CCBr') == 0 and has_ez('BrCC(C)=CCBr')
tot = sum(nstereo(s) for s in iso); print('立体異性体の合計', tot); assert tot == 6
# (3) 2 分子目の付加
from rdkit.Chem import AllChem
r12 = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1](Br)[C:2]Br')
second = set()
for s in iso:
    for ps in r12.RunReactants((Chem.MolFromSmiles(s),)):
        p = ps[0]; Chem.SanitizeMol(p); second.add(Chem.MolToSmiles(p))
print('2 分子目', second); assert second == {canon('BrCC(C)(Br)C(Br)CBr')}
assert nstereo('BrCC(C)(Br)C(Br)CBr') == 4
# (4) ペンタジエン
a13 = adducts('C=CC=CC'); a14 = adducts('C=CCC=C')
print('1,3-ペンタジエン', a13, '1,4-ペンタジエン', a14)
assert len(a13) == 3 and len(a14) == 1
print('すべて OK')
