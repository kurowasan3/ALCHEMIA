# 検証スクリプト THE-BND-093   python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-BND-093/item.json'
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
# 5CB と 2,2-ジメチルプロパンの分子式
assert CalcMolFormula(Chem.MolFromSmiles('CCCCCc1ccc(-c2ccc(C#N)cc2)cc1')) == 'C18H19N'
assert CalcMolFormula(Chem.MolFromSmiles('CC(C)(C)C')) == 'C5H12'
# 分類表：(位置, 向き) -> 状態，流動性は位置の秩序がないもの
table = {'ア': (1, 1, '結晶'), 'イ': (1, 0, '柔粘性結晶'), 'ウ': (0, 1, '液晶'), 'エ': (0, 0, '液体')}
fluid = [k for k, (p, o, _) in table.items() if p == 0]
assert fluid == ['ウ', 'エ']
print('分類 OK')

print('検証OK')
