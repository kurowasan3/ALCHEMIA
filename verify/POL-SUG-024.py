# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/out/items/POL-SUG-024/item.json'
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
        out.append(s[i + 4:j - 1]); i = j


eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
# 小問数と答えの番号の一致
nq = it['problem'].count('\\item')
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問', nq, '答え', na)


def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import sympy as sp
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
Cu2O = 2 * 64 + 16
assert Cu2O == 144
nmol = 0.288 / Cu2O
assert sig(nmol, 3) == 0.00200
Mav = 3.28 / nmol
n = (Mav - 18) / 162
print('平均分子量', Mav, '重合度', n)
assert round(n) == 10 and abs(n - 10) < 0.05
g = 10 * 0.00200
assert sig(g * Cu2O, 3) == 2.88
ratio = (162 * 10 + 18) / (162 * 1000 + 18)
print('比', ratio); assert sig(ratio, 1) == 0.01
# グルコース 10 単位の直鎖（α-1,4）の分子式
# グルコース 10 単位の直鎖（式だけの確認なので立体は省略）
chain = 'OCC1OC(O)C(O)C(O)C1' + 'OC1OC(CO)C(O)C(O)C1' * 8 + 'OC1OC(CO)C(O)C(O)C1O'
mm = Chem.MolFromSmiles(chain)
print(CalcMolFormula(mm), MolWt(mm))
assert CalcMolFormula(mm) == 'C60H102O51'
assert 60 * 12 + 102 * 1.0 + 51 * 16 == 162 * 10 + 18
print('POL-SUG-024 検証OK')
