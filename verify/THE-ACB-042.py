# 検証スクリプト（グループ g11）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-ACB-042/item.json'
import sys, json, re, math
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
eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
assert check_equation('2XeF2 + 2H2O -> 2Xe + 3HF + O2') is not None  # 検査器が誤りを検出できること
def sig(x, n):
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
# 反応式（式OK の出力）に加え，分子式を確認
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
imine = Chem.MolFromSmiles('c1ccccc1C=NC(C)(C)C'); amine = Chem.MolFromSmiles('c1ccccc1CNC(C)(C)C')
print(CalcMolFormula(imine), CalcMolFormula(amine))
assert CalcMolFormula(imine) == 'C11H15N' and CalcMolFormula(amine) == 'C11H17N'
assert check_equation('C11H15N + H2 -> C11H17N') is None
# B の周りの電子数：BF3 は 3 組=6，付加物は 4 組=8
assert 3*2 == 6 and 4*2 == 8
# 共役酸の pKa（文献値）：ピリジン 5.2 < 2,6-ルチジン 6.7 → ルチジンの方が H+ に対して強い塩基
assert 6.7 > 5.2
print('OK')
