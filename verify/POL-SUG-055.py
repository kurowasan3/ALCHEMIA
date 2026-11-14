# 検証スクリプト POL-SUG-055  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/POL-SUG-055/item.json'
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|\$', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式の検査', n, '本 OK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

import numpy as np, itertools
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
POS = {'U': np.array([0, 1, -1.]), 'D': np.array([0, -1, -1.]), 'L': np.array([-1, 0, 1.]), 'R': np.array([1, 0, 1.])}
S = -1  # 較正済み（下で確認）
def center_smiles(groups):
    """groups: dict 位置->基 ('CHO','OH','CH2OH','H') の 1 炭素の Fischer 式 → SMILES"""
    frag = {'CHO': 'C=O', 'OH': 'O', 'CH2OH': 'CO'}
    pos = {g: p for p, g in groups.items()}
    order = ['CHO', 'H', 'OH', 'CH2OH']   # SMILES: [C@?H](C=O)... 書き方を固定
    # SMILES "O=C[C@?H](O)CO" の近接順は CHO, H, OH, CH2OH
    v = np.linalg.det(np.array([POS[pos['H']] - POS[pos['CHO']], POS[pos['OH']] - POS[pos['CHO']], POS[pos['CH2OH']] - POS[pos['CHO']]]))
    tg = '@' if v*S > 0 else '@@'
    return Chem.CanonSmiles(f'O=C[C{tg}H](O)CO')
def chain_smiles(sides):
    out = 'O=C'
    for x in sides:
        h = 'L' if x == 'R' else 'R'
        v = np.linalg.det(np.array([POS[h] - POS['U'], POS[x] - POS['U'], POS['D'] - POS['U']]))
        out += f"[C{'@' if v*S > 0 else '@@'}H](O)"
    return Chem.CanonSmiles(out + 'CO')
def cip(sm):
    m = Chem.MolFromSmiles(sm); rdCIPLabeler.AssignCIPLabels(m)
    return ''.join(a.GetProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode'))
# 較正
Dg = center_smiles({'U': 'CHO', 'R': 'OH', 'D': 'CH2OH', 'L': 'H'})
assert cip(Dg) == 'R' and chain_smiles('R') == Dg
assert cip(chain_smiles('RLRR')) == 'RSRR'
Lg = center_smiles({'U': 'CHO', 'R': 'H', 'D': 'CH2OH', 'L': 'OH'})
a = center_smiles({'U': 'CH2OH', 'R': 'H', 'D': 'CHO', 'L': 'OH'})
b = center_smiles({'U': 'CHO', 'R': 'H', 'D': 'CH2OH', 'L': 'OH'})
c = center_smiles({'U': 'H', 'R': 'CHO', 'D': 'OH', 'L': 'CH2OH'})
print('a', a == Dg, 'b', b == Lg, 'c', c == Lg)
assert a == Dg and b == Lg and c == Lg
names = {'RRR': 'D-ribose', 'LRR': 'D-arabinose', 'RLR': 'D-xylose', 'LLR': 'D-lyxose'}
known = {'D-ribose': 'RRR', 'D-arabinose': 'SRR', 'D-xylose': 'RSR', 'D-lyxose': 'SSR'}
smi = {}
for sd, nm in names.items():
    smi[nm] = chain_smiles(sd); print(nm, cip(smi[nm])); assert cip(smi[nm]) == known[nm]
def mirror(sd): return ''.join('L' if x == 'R' else 'R' for x in sd)
for sd, nm in list(names.items()):
    smi['L-' + nm[2:]] = chain_smiles(mirror(sd))
allst = {chain_smiles(''.join(t)) for t in itertools.product('LR', repeat=3)}
assert len(allst) == 8
xyl = 'RLR'
epi = []
for i in range(3):
    s = list(xyl); s[i] = 'L' if s[i] == 'R' else 'R'
    sm = chain_smiles(''.join(s)); epi.append([k for k, v in smi.items() if v == sm][0])
print('epimers', epi); assert epi == ['D-lyxose', 'D-ribose', 'L-arabinose']
dia = allst - {smi['D-xylose'], smi['L-xylose']}; assert len(dia) == 6

print('ALL OK')
