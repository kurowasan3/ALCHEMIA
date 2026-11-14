# 検証スクリプト POL-RES-016  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/POL-RES-016/item.json'
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
from rdkit.Chem.Descriptors import MolWt
W = {'C': 12, 'H': 1.0, 'O': 16}
def M(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); return sum(W[a.GetSymbol()] for a in m.GetAtoms())
assert CalcMolFormula(Chem.MolFromSmiles('COC(=O)c1ccc(C(=O)OC)cc1')) == 'C10H10O4' and M('COC(=O)c1ccc(C(=O)OC)cc1') == 194
assert CalcMolFormula(Chem.MolFromSmiles('OCCOC(=O)c1ccc(C(=O)OCCO)cc1')) == 'C12H14O6' and M('OCCOC(=O)c1ccc(C(=O)OCCO)cc1') == 254
assert 10*12 + 8*1.0 + 4*16 == 192 and M('OCCO') == 62 and M('CO') == 32
n = 96.0/192
assert (round(n*194, 1), round(n*62, 1), round(2*n*32, 1), round(n*254)) == (97.0, 31.0, 32.0, 127)
assert 96.0 + 32.0 == 97.0 + 31.0
import json
P = json.dumps(json.load(open('/home/claude/kurowasan3/alchemia/content/chem/polymers.json', encoding='utf-8')), ensure_ascii=False)
assert '熱硬化' in P
print('PET 解重合の量的関係 OK')

print('ALL OK')
