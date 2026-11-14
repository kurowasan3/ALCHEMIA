# ORG-DET-018 検証：CH3-X, C6H5-X, H-X の組と性質（tests.json の判定を使う）
import json
T = {t['id']: t for t in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/tests.json'))['tests']}
def pos(test, cid):
    t = T[test]
    if cid in t['pos']: return True
    if cid in t['scope']: return False
    return None
# 図鑑にない物質は高校教科書の標準的知識で補う
EXTRA = {
 ('sodium', 'water'): True, ('silver-mirror', 'water'): False, ('fecl3', 'water'): False, ('nahco3', 'water'): False,
 ('sodium', 'methylamine'): False, ('sodium', 'ammonia'): False, ('sodium', 'chlorobenzene'): False, ('sodium', 'hcl'): None,
 ('silver-mirror', 'methylamine'): False, ('silver-mirror', 'ammonia'): False, ('silver-mirror', 'chloromethane'): False,
 ('silver-mirror', 'chlorobenzene'): False, ('silver-mirror', 'hcl'): False, ('silver-mirror', 'propene'): False,
 ('silver-mirror', 'styrene'): False, ('silver-mirror', 'methane'): False,
}
X = {
 'ア-OH': ['methanol', 'phenol', 'water'],
 'イ-CHO': ['acetaldehyde', 'benzaldehyde', 'formaldehyde'],
 'ウ-COOH': ['acetic-acid', 'benzoic-acid', 'formic-acid'],
 'エ-NH2': ['methylamine', 'aniline', 'ammonia'],
 'オ-Cl': ['chloromethane', 'chlorobenzene', 'hcl'],
 'カ-CH=CH2': ['propene', 'styrene', 'ethylene'],
 'キ-CH3': ['ethane', 'toluene', 'methane'],
}
def p(test, c):
    v = pos(test, c) if test in T else None
    if v is None: v = EXTRA.get((test, c))
    return v
# (1) 3つとも Na と反応
r1 = [k for k, cs in X.items() if all(p('sodium', c) for c in cs)]
# (2) 3つとも銀鏡 / (3) 1つだけ銀鏡
sm = {k: [p('silver-mirror', c) for c in cs] for k, cs in X.items()}
r2 = [k for k, v in sm.items() if all(v)]; r3 = [k for k, v in sm.items() if sum(bool(x) for x in v) == 1]
# (4) 3つとも NaHCO3 で CO2
r4 = [k for k, cs in X.items() if all(p('nahco3', c) for c in cs)]
# (5) 1つだけ FeCl3 呈色
r5 = [k for k, cs in X.items() if sum(bool(p('fecl3', c)) for c in cs) == 1]
print(r1, r2, r3, r4, r5)
assert r1 == ['ア-OH', 'ウ-COOH'] and r2 == ['イ-CHO'] and r3 == ['ウ-COOH'] and r4 == ['ウ-COOH'] and r5 == ['ア-OH']
# (6) 塩酸と塩：アミン類・アンモニア（教科書知識）／(7) 付加重合：C=C をもつ単量体／(8) 水に溶けて強酸：HCl のみ
from rdkit import Chem
smi = {'methylamine': 'CN', 'aniline': 'Nc1ccccc1', 'ammonia': 'N', 'propene': 'C=CC', 'styrene': 'C=Cc1ccccc1', 'ethylene': 'C=C'}
assert all(Chem.MolFromSmiles(smi[c]).HasSubstructMatch(Chem.MolFromSmarts('[NX3;H2,H3;!$(NC=O)]')) for c in X['エ-NH2'])
assert all(Chem.MolFromSmiles(smi[c]).HasSubstructMatch(Chem.MolFromSmarts('C=C')) for c in X['カ-CH=CH2'])
print('ORG-DET-018 OK')
