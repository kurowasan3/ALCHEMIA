# 検証スクリプト ORG-ARO-049  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ARO-049/item.json'
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
def nC(smi):
    m = Chem.MolFromSmiles(smi); r = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    return len({r[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol() == 'C'})
exp = {'Cc1ccccc1': 5, 'Cc1ccccc1O': 7, 'Cc1cccc(O)c1': 7, 'Cc1ccc(O)cc1': 5, 'OCc1ccccc1': 5, 'COc1ccccc1': 5,
       'OC(=O)c1ccccc1C(=O)O': 4, 'OC(=O)c1cccc(C(=O)O)c1': 5, 'OC(=O)c1ccc(C(=O)O)cc1': 3, 'O=C1OC(=O)c2ccccc21': 4}
for s, n in exp.items(): assert nC(s) == n, (s, nC(s), n)
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C7H8O' for s in list(exp)[1:6])
print('炭素の種類の数 OK')

print('ALL OK')
