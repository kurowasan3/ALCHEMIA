# POL-PRT-046 検証：トリプトファン→メラトニンの分子式・不斉炭素・溶解性
# --- 共通：item.json の \ce{} 反応式のつり合い検査（単独で動く） ---
import json, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
def _ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k, d = j + 4, 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k
def check_item_equations(path):
    it = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for ce in _ce_blocks(it[key]):
            arrows = re.findall(r'<=>|->|<-', ce)
            if len(arrows) != 1 or '\\boxed' in ce or '\\text' in ce or '\\mathrm' in ce:
                continue
            r = check_equation(ce)
            assert r is None, (key, ce, r)
            n += 1
    print('反応式', n, '本 OK')
    return it
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import Counter
check_item_equations('/home/claude/chemx/todai/out/items/POL-PRT-046/item.json')
S = {'Trp': 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O', 'X': 'N[C@@H](Cc1c[nH]c2ccc(O)cc12)C(=O)O',
     'Y': 'NCCc1c[nH]c2ccc(O)cc12', 'Z': 'CC(=O)NCCc1c[nH]c2ccc(O)cc12', 'Mel': 'CC(=O)NCCc1c[nH]c2ccc(OC)cc12'}
F = {'Trp': 'C11H12N2O2', 'X': 'C11H12N2O3', 'Y': 'C10H12N2O', 'Z': 'C12H14N2O2', 'Mel': 'C13H16N2O2'}
M = {k: Chem.MolFromSmiles(v) for k, v in S.items()}
for k in S: assert CalcMolFormula(M[k]) == F[k], k
def comp(m): m = Chem.AddHs(m); return Counter(a.GetSymbol() for a in m.GetAtoms())
seq = ['Trp', 'X', 'Y', 'Z', 'Mel']
for a, b in zip(seq, seq[1:]):
    d = comp(M[b]); d.subtract(comp(M[a])); print(a, '->', b, {k: v for k, v in d.items() if v})
chiral = [k for k in seq if Chem.FindMolChiralCenters(M[k], includeUnassigned=True, useLegacyImplementation=False)]
print('chiral', chiral); assert chiral == ['Trp', 'X']
pat = {'amine': '[NX3;H2][CX4]', 'cooh': 'C(=O)[OH]', 'phenol': 'c[OH]'}
both = []
for k in seq:
    h = {p: M[k].HasSubstructMatch(Chem.MolFromSmarts(s)) for p, s in pat.items()}
    if h['amine'] and (h['cooh'] or h['phenol']): both.append(k)
    print(k, h)
assert both == ['Trp', 'X', 'Y']
print('OK')
