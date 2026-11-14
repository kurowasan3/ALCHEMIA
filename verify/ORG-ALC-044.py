# 検証スクリプト ORG-ALC-044   python3 -I この.py
ITEM = '/home/claude/chemk/out/items/ORG-ALC-044/item.json'
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
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Descriptors
def f(s): return CalcMolFormula(Chem.MolFromSmiles(s))
# (2) プロピル + CO2 → ブタン酸（C3 + 1C）
assert f('CCCC(=O)O') == 'C4H8O2'
# 2-ブロモプロパン由来の R=(CH3)2CH にCOOH → 2-メチルプロパン酸
assert Chem.CanonSmiles('CC(C)C(=O)O') == Chem.CanonSmiles('OC(=O)C(C)C')
# (3) D 入りの式を H に置き換えてつり合いを確認（重水素は chemlib で読めないため）
e = 'CH3CH(CH3)MgBr + D2O -> CH3CHDCH3 + Mg(OD)Br'.replace('D', 'H')
assert check_equation(e) is None, check_equation(e)
print('D 置換式 OK')
# (4) フェニル + エチレンオキシド → 2-フェニルエタノール（C6→C8）
assert f('OCCc1ccccc1') == 'C8H10O'
# (5) メチルオキシランの開環：CH2 側攻撃 → 2-ブタノール，CH 側攻撃 → 2-メチル-1-プロパノール
a = Chem.CanonSmiles('CC(O)CC'); b = Chem.CanonSmiles('CC(C)CO')
assert f(a) == f(b) == 'C4H10O'
# ヨードホルム陽性構造 CH3CH(OH)- の有無
patt = Chem.MolFromSmarts('[CH3][CH1]([OH])')
assert Chem.MolFromSmiles(a).HasSubstructMatch(patt) and not Chem.MolFromSmiles(b).HasSubstructMatch(patt)
print('ヨードホルム判定 OK')

print('検証OK')
