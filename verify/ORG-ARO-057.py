# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ARO-057/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
print('小問の数', np_)

from math import log10, floor, sqrt
def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Descriptors
smi = {'ベンジルアルコール':'OCc1ccccc1','o-クレゾール':'Cc1ccccc1O','m-クレゾール':'Cc1cccc(O)c1','p-クレゾール':'Cc1ccc(O)cc1','アニソール':'COc1ccccc1'}
for n, s in smi.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C7H8O', n
    print(n, CalcMolFormula(m), round(Descriptors.MolWt(m),2))
# 芳香族 C7H8O（ベンゼン環1つ）を列挙：フェノール類は 3 種類
phen = set(Chem.CanonSmiles(s) for s in ['Cc1ccccc1O','Cc1cccc(O)c1','Cc1ccc(O)cc1','Oc1ccccc1C','Oc1cc(C)ccc1'])
assert len(phen) == 3
M = 7*12 + 8*1.0 + 16
assert abs(M - 108) < 1e-9
ntot = 5.40 / M
nH2 = 0.336 / 22.4
nOH = 2 * nH2
nY = 1.08 / M
nX = nOH - nY; nZ = ntot - nOH
print('ntot', ntot, 'nX', nX, 'nY', nY, 'nZ', nZ)
assert sig(nX, 3) == 0.0200 and sig(nY, 3) == 0.0100 and sig(nZ, 3) == 0.0200
for n, p in [(nX, 40.0), (nY, 20.0), (nZ, 40.0)]:
    assert sig(n/ntot*100, 3) == p
V5 = nX / 2 * 22.4 * 1000
print('(5) H2 mL', V5); assert sig(V5, 3) == 224
# 酸の強さ：クレゾール pKa 約 10.3 ＞ 炭酸 pKa1 6.35 → NaHCO3 とは反応しない
K = 10**-10.3 / 10**-6.35; print('クレゾール + HCO3- の K', K); assert K < 1e-3
print('ORG-ARO-057 OK')
