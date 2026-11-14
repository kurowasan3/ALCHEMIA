# 検証スクリプト POL-PRT-035  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/POL-PRT-035/item.json'
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
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from fractions import Fraction as F
AW = {'H':F(1),'C':F(12),'N':F(14),'O':F(16),'S':F(32)}
def mw(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    return sum(AW[a.GetSymbol()] for a in m.GetAtoms())
P = 'NCC(=O)NC(CS)C(=O)O'; Q = 'CC(N)C(=O)NC(CS)C(=O)O'
assert mw(P) == 178 and mw(Q) == 192
PP = 'NCC(=O)NC(CSSCC(NC(=O)CN)C(=O)O)C(=O)O'
PQ = 'NCC(=O)NC(CSSCC(NC(=O)C(C)N)C(=O)O)C(=O)O'
QQ = 'CC(N)C(=O)NC(CSSCC(NC(=O)C(C)N)C(=O)O)C(=O)O'
assert (mw(PP), mw(PQ), mw(QQ)) == (354, 368, 382)
assert CalcMolFormula(Chem.MolFromSmiles('NC(CSSCC(N)C(=O)O)C(=O)O')) == 'C6H12N2O4S2'
n = lambda k: k*(k+1)//2
assert n(2) == 3 and n(3) == 6
mass = F(178,100)/178/2*354
assert mass == F(177,100)
print('分子量・組合せ・質量 OK')

print('ALL OK')
