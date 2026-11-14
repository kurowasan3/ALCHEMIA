# 検証スクリプト（グループ g12）  python3 -I このファイル
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

def check_all_eqs(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)

def check_parts(it):
    n = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n + 1)], (n, nums)
    print('小問数OK:', n)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

ITEM = '/home/claude/chemk/out/items/ORG-ARO-096/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_all_eqs(it); check_parts(it)
# 追加の式（解説外の確認）
for e in ['C6H5COC2H5 + NH2OH -> C6H5C(=NOH)C2H5 + H2O',
          'C6H5NHCOC2H5 + H2O + HCl -> C6H5NH3Cl + C2H5COOH',
          'C5H8O + NH2OH -> C5H9NO + H2O']:
    assert check_equation(e) is None, e
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
def f(s): return F(Chem.MolFromSmiles(s))
P = 'CC(=O)Nc1ccccc1'; Q = 'CNC(=O)c1ccccc1'
oxX = r'C/C(=N\O)c1ccccc1'
assert f(P) == f(Q) == 'C8H9NO'
assert f(oxX) == 'C8H9NO'          # オキシムとアミドは異性体（転位で原子の出入りなし）
ox5 = 'ON=C1CCCC1'; lac = 'O=C1CCCCN1'
assert f(ox5) == f(lac) == 'C5H9NO'
mlac = Chem.MolFromSmiles(lac); ri = mlac.GetRingInfo().AtomRings()
assert len(ri) == 1 and len(ri[0]) == 6
# (5) N-フェニルプロパンアミド + H2O -> アニリン + プロピオン酸（原子の保存）
from collections import Counter
def cnt(s):
    m = Chem.AddHs(Chem.MolFromSmiles(s)); return Counter(a.GetSymbol() for a in m.GetAtoms())
assert cnt('CCC(=O)Nc1ccccc1') + cnt('O') == cnt('Nc1ccccc1') + cnt('CCC(=O)O')
assert f('CCC(=O)c1ccccc1') == 'C9H10O'
print('RDKit OK')
print('ALL OK')
