# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/outh/items/POL-PRT-050/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 小問の数と答えの数の一致
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 0), (nprob, nans)
print('小問数', nprob, '答え', nans)

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import permutations
cyc = Chem.MolFromSmiles('O=C1CNC(=O)CNC(=O)C(C)NC(=O)C(CCCCN)N1')
lin = Chem.MolFromSmiles('NC(C)C(=O)NCC(=O)NCC(=O)NC(CCCCN)C(=O)O')
assert CalcMolFormula(cyc) == 'C13H23N5O4' and CalcMolFormula(lin) == 'C13H25N5O5'
Mres = {'G': 75-18, 'A': 89-18, 'K': 146-18}
assert 2*Mres['G'] + Mres['A'] + Mres['K'] == 313
water = (15.0 + 8.9 + 14.6 - 31.3)/18; res = 15.0/75 + 8.9/89 + 14.6/146
assert abs(water - res) < 1e-9 and abs(water - 0.4) < 1e-9
cycles = set()
for p in set(permutations('GGAK')):
    rots = [p[i:] + p[:i] for i in range(4)]
    cycles.add(min(rots))
print(cycles); assert len(cycles) == 3
ok = []
for c in cycles:
    i = c.index('K'); s = c[i:] + c[:i]          # K から向きにそって読む
    pairs = [s[j] + s[(j+1) % 4] for j in range(4)]
    product = s[1:] + s[:1]                      # K の C 側で開環：K の次から始まり K で終わる
    if 'GG' in pairs and product[0] == 'A': ok.append(product)
print(ok); assert ok == [('A', 'G', 'G', 'K')]

print('ALL OK')
