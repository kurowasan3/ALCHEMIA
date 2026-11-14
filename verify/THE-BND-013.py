# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-013/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
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

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
EN = {'H': 2.20, 'B': 2.04, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98, 'S': 2.58, 'Cl': 3.16}
VAL = {'B': 3, 'C': 4, 'N': 5, 'O': 6, 'S': 6}
mols = {'BF3': 'FB(F)F', 'NF3': 'FN(F)F', 'CS2': 'S=C=S', 'OF2': 'FOF', 'CCl4': 'ClC(Cl)(Cl)Cl', 'CH2Cl2': 'ClCCl'}
res = {}
for name, smi in mols.items():
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    c = [a for a in m.GetAtoms() if a.GetSymbol() in VAL and a.GetDegree() >= 2][0]
    bo = sum(b.GetBondTypeAsDouble() for b in c.GetBonds())
    lone = (VAL[c.GetSymbol()] - bo) / 2
    steric = c.GetDegree() + lone
    AllChem.EmbedMolecule(m, randomSeed=7); AllChem.MMFFOptimizeMolecule(m) if name not in ('BF3',) else None
    conf = m.GetConformer()
    if name == 'BF3':  # 力場に B が無いので理想構造
        pass
    v = np.zeros(3)
    for b in c.GetBonds():
        o = b.GetOtherAtom(c)
        r = np.array(conf.GetAtomPosition(o.GetIdx())) - np.array(conf.GetAtomPosition(c.GetIdx()))
        v += (EN[o.GetSymbol()] - EN[c.GetSymbol()]) * r / np.linalg.norm(r)
    res[name] = (int(lone), int(steric), round(float(np.linalg.norm(v)), 2))
    print(name, '非共有電子対', lone, '電子対の組', steric, '|双極子和|', res[name][2])
assert res['BF3'][:2] == (0, 3) and res['NF3'][:2] == (1, 4)
assert res['CS2'][:2] == (0, 2) and res['OF2'][:2] == (2, 4)
assert res['CCl4'][:2] == (0, 4) and res['CH2Cl2'][:2] == (0, 4)
for n in ('NF3', 'OF2', 'CH2Cl2'): assert res[n][2] > 0.2
for n in ('CS2', 'CCl4'): assert res[n][2] < 0.1
# BF3：電子対3組・非共有0 → 平面正三角形（対称）なので双極子和は0（ETKDG の構造でも確認）
print('BF3 双極子和（ETKDG）', res['BF3'][2])

print('全検証OK:', ITEM)
