# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/ORG-CBX-060/item.json'
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

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ['src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer', 'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig']:
    assert k in it, k
# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
def iw(s):
    m = Chem.AddHs(Chem.MolFromSmiles(s)); w = {1: 1, 6: 12, 8: 16}
    return sum(w[a.GetAtomicNum()] for a in m.GetAtoms())
G, Su, Hb = 'OCCO', 'OC(=O)CCC(=O)O', 'OCCCC(=O)O'
J = 'OCCOC(=O)CCC(=O)O'; K = 'O=C1CCC(=O)OCCO1'; L = 'OC(=O)CCC(=O)OCCOC(=O)CCC(=O)O'; M = 'O=C1CCCO1'
assert (iw(G), iw(Su), iw(Hb)) == (62, 118, 104)
assert [iw(x) for x in (J, K, L, M)] == [162, 144, 262, 86]
water = [(iw(G)+iw(Su)-iw(J))/18, (iw(G)+iw(Su)-iw(K))/18, (iw(G)+2*iw(Su)-iw(L))/18, (iw(Hb)-iw(M))/18]
print(water); assert water == [1, 2, 2, 1]
mols = [Chem.MolFromSmiles(x) for x in (J, K, L, M)]
ri = [mm.GetRingInfo().AtomRings() for mm in mols]
print([[len(r) for r in x] for x in ri]); assert [len(r) for r in ri[1]] == [8] and [len(r) for r in ri[3]] == [5] and not ri[0] and not ri[2]
cooh = Chem.MolFromSmarts('C(=O)[OX2H1]')
print([len(Chem.MolFromSmiles(x).GetSubstructMatches(cooh)) for x in (J, K, L, M)])
assert [len(Chem.MolFromSmiles(x).GetSubstructMatches(cooh)) for x in (J, K, L, M)] == [1, 0, 2, 0]
assert round(2.62/262*4*40, 2) == 1.60
print('OK')
