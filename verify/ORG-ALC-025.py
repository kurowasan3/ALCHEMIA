# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ALC-025/item.json'
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
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

# ---- 計算
from sympy import Rational as R, symbols, solve

# ---- RDKit：官能基による検出反応の判定と水素付加の生成物
from rdkit import Chem
cpd = {'A': 'C=CCO', 'B': 'C=CC=O', 'C': 'C=CC(=O)O', 'D': 'CCCO', 'E': 'CCC=O'}
pat = {'a': '[CX3]=[CX3]', 'b': '[CX3H1](=O)[#6]', 'c': 'C(=O)[OH1]'}
res = {k: sorted(n for n, s in cpd.items() if Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts(p))) for k, p in pat.items()}
print(res)
assert res == {'a': ['A', 'B', 'C'], 'b': ['B', 'E'], 'c': ['C']}
# Ni 触媒の水素付加：C=C とアルデヒドの C=O を還元，COOH は還元しない（高校の標準的な扱い）
def hydrogenate(s):
    rw = Chem.RWMol(Chem.MolFromSmiles(s)); nH2 = 0
    for b in rw.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if b.GetBondTypeAsDouble() != 2: continue
        if a1.GetSymbol() == a2.GetSymbol() == 'C':
            b.SetBondType(Chem.BondType.SINGLE); nH2 += 1
        else:
            c = a1 if a1.GetSymbol() == 'C' else a2
            if not any(x.GetSymbol() == 'O' and rw.GetBondBetweenAtoms(c.GetIdx(), x.GetIdx()).GetBondTypeAsDouble() == 1 for x in c.GetNeighbors()):
                b.SetBondType(Chem.BondType.SINGLE); nH2 += 1
    m = rw.GetMol(); Chem.SanitizeMol(m); return Chem.MolToSmiles(m), nH2
assert hydrogenate('C=CC=O') == (Chem.CanonSmiles('CCCO'), 2)
assert hydrogenate('C=CC(=O)O') == (Chem.CanonSmiles('CCC(=O)O'), 1)
assert hydrogenate('C=CCO') == (Chem.CanonSmiles('CCCO'), 1)
print('OK')
