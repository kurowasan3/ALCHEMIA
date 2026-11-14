# 検証スクリプト（グループ g03）  python3 -I この.py

import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_all_ce(it):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    print('検査した反応式:', n)

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
it = json.load(open('/home/claude/chem5/out/items/ORG-CBX-030/item.json', encoding='utf-8'))
check_all_ce(it)
# C4H8O2 の鎖式カルボン酸・エステルの列挙：RCOOR' と RCOOH（R は H またはアルキル基）
alkyl = {0: ['[H]'], 1: ['C'], 2: ['CC'], 3: ['CCC', 'C(C)C']}
acids, esters = set(), set()
for r in alkyl[3]:
    acids.add(Chem.MolToSmiles(Chem.MolFromSmiles(r + 'C(=O)O') if r != 'C(C)C' else Chem.MolFromSmiles('CC(C)C(=O)O')))
for nr in range(0, 3):
    for r in alkyl[nr]:
        for rp in alkyl[3 - nr]:
            if rp == '[H]': continue
            head = 'O=C' if r == '[H]' else r + 'C(=O)'
            s = (head + 'O' + rp) if r != '[H]' else ('O=CO' + rp)
            if rp == 'C(C)C': s = (head + 'OC(C)C') if r != '[H]' else 'O=COC(C)C'
            esters.add(Chem.MolToSmiles(Chem.MolFromSmiles(s)))
print('カルボン酸', sorted(acids)); print('エステル', sorted(esters))
assert len(acids) == 2 and len(esters) == 4
for s in acids | esters:
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C4H8O2'; assert round(Descriptors.MolWt(m)) == 88
# 示性式との対応
ans = {'CCCC(=O)O', 'CC(C)C(=O)O', 'CCCOC=O', 'CC(C)OC=O', 'CCOC(C)=O', 'CCC(=O)OC'}
assert {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in ans} == acids | esters
formates = [s for s in esters if Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CH1](=O)O[#6]'))]
assert len(formates) == 2
# 沸点（文献値，℃）
bp = {'CCCC(=O)O': 164, 'CC(C)C(=O)O': 155, 'CCCOC=O': 81, 'CC(C)OC=O': 68, 'CCOC(C)=O': 77, 'CCC(=O)OC': 80}
bp = {Chem.MolToSmiles(Chem.MolFromSmiles(k)): v for k, v in bp.items()}
assert min(bp[a] for a in acids) > max(bp[e] for e in esters)
assert sorted(bp.values()) == [68, 77, 80, 81, 155, 164]
print('すべて一致')
