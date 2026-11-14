"""ORG-DET-027 検証：異性体の組の分子式，各操作の陽性・陰性表，反応式"""
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
pairs = {
 'ア': ('CCO', 'COC'),
 'イ': ('CCC=O', 'CC(C)=O'),
 'ウ': ('CCCO', 'CC(C)O'),
 'エ': ('OC(=O)/C=C\\C(=O)O', 'OC(=O)/C=C/C(=O)O'),
 'オ': ('Cc1ccccc1O', 'OCc1ccccc1'),
 'カ': ('CC(=O)O', 'O=COC'),
}
for k, (a, b) in pairs.items():
    fa, fb = CalcMolFormula(Chem.MolFromSmiles(a)), CalcMolFormula(Chem.MolFromSmiles(b))
    assert fa == fb, (k, fa, fb); print(k, fa)
# 陽性表（1=陽性）: A Na, B NaHCO3, C 銀鏡, D ヨードホルム, E FeCl3, F 臭素水(C=C付加・フェノール置換のみ), G 加熱で酸無水物
T = {
 'エタノール': 'A D', 'ジメチルエーテル': '', 'プロパナール': 'C', 'アセトン': 'D',
 '1-プロパノール': 'A', '2-プロパノール': 'A D', 'マレイン酸': 'A B F G', 'フマル酸': 'A B F',
 'o-クレゾール': 'A E F', 'ベンジルアルコール': 'A', '酢酸': 'A B', 'ギ酸メチル': 'C',
}
P = {'ア': ('エタノール', 'ジメチルエーテル'), 'イ': ('プロパナール', 'アセトン'), 'ウ': ('1-プロパノール', '2-プロパノール'),
     'エ': ('マレイン酸', 'フマル酸'), 'オ': ('o-クレゾール', 'ベンジルアルコール'), 'カ': ('酢酸', 'ギ酸メチル')}
expect = {'ア': 'AD', 'イ': 'CD', 'ウ': 'D', 'エ': 'G', 'オ': 'EF', 'カ': 'ABC'}
for k, (a, b) in P.items():
    usable = ''.join(sorted(set(T[a].split()) ^ set(T[b].split())))
    print(k, usable); assert usable == expect[k]
# tests.json との照合（記載のあるもの）
tj = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/tests.json'))['tests']
pos = {t['id']: set(t['pos']) for t in tj}
assert 'ethanol' in pos['iodoform'] and 'propanol-2' in pos['iodoform'] and 'acetone' in pos['iodoform']
assert 'acetic-acid' not in pos['iodoform'] and 'propanol-1' not in pos['iodoform']
assert 'o-cresol' in pos['fecl3'] and 'benzyl-alcohol' not in pos['fecl3']
assert 'o-cresol' in pos['bromine'] and 'acetic-acid' in pos['nahco3'] and 'benzyl-alcohol' in pos['sodium']
for eq in ['C4H4O4 -> C4H2O3 + H2O', '2C3H7OH + 2Na -> 2C3H7ONa + H2']:
    r = check_equation(eq); print(eq, r or 'OK'); assert r is None
# item 内の ce 式をすべて検査
it = json.load(open('/home/claude/chem5/out/items/ORG-DET-027/item.json'))
for key in ('problem', 'answer', 'solution'):
    for m in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[key]):
        if re.search(r'->|<=>', m) and '\\text' not in m and '\\mathrm' not in m:
            r = check_equation(m); assert r is None, (m, r)
print('ALL OK')
