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

ITEM = '/home/claude/chemk/outh/items/ORG-HYD-061/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_all_eqs(it); check_parts(it)
for e in ['C5H11CH=CHCHO + (C6H5)3P=CHCOOC2H5 -> C5H11CH=CHCH=CHCOOC2H5 + (C6H5)3PO',
          'OHCCH=CHCOOC2H5 + (C6H5)3P=CHC5H11 -> C5H11CH=CHCH=CHCOOC2H5 + (C6H5)3PO']:
    assert check_equation(e) is None, e
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
T = r'CCOC(=O)/C=C/C=C\CCCCC'
mT = Chem.MolFromSmiles(T)
assert F(mT) == 'C12H20O2'
st = {}
for b in mT.GetBonds():
    if b.GetStereo() != Chem.BondStereo.STEREONONE:
        st[(b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx())] = str(b.GetStereo())
print('T の二重結合の立体:', st)
flat = Chem.MolFromSmiles('CCOC(=O)C=CC=CCCCCC')
isos = {Chem.MolToSmiles(m) for m in EnumerateStereoisomers(flat)}
assert len(isos) == 4; print('立体異性体の数 =', len(isos))

# モデルによる全組み合わせの判定
# アルデヒド: (名前, 骨格側の位置, もとの C=C の形)  骨格: C1(エステル)…C10
ald = {'ヘキサナール': ('C5', None), 'トランス-2-オクテナール': ('C3', ('C4=C5', 'E')),
       'シス-2-オクテナール': ('C3', ('C4=C5', 'Z')), 'トランス-4-オキソ-2-ブテン酸エチル': ('C4', ('C2=C3', 'E'))}
# 臭化物: (イリド炭素の位置, 安定か, もとの C=C)
brm = {'ブロモ酢酸エチル': ('C2', True, None), '1-ブロモヘキサン': ('C5', False, None),
       'トランス-4-ブロモ-2-ブテン酸エチル': ('C4', True, ('C2=C3', 'E')),
       'シス-4-ブロモ-2-ブテン酸エチル': ('C4', True, ('C2=C3', 'Z'))}
# 骨格が成り立つ組み合わせ：アルデヒド炭素とイリド炭素が隣り合い，両方合わせて C1〜C10 になる
valid_pairs = {('C3', 'C2'): 'C2=C3', ('C4', 'C5'): 'C4=C5', ('C5', 'C4'): 'C4=C5'}
# 各原料の炭素範囲で骨格の重なり・不足が無いか
span_ald = {'ヘキサナール': (5, 10), 'トランス-2-オクテナール': (3, 10), 'シス-2-オクテナール': (3, 10), 'トランス-4-オキソ-2-ブテン酸エチル': (1, 4)}
span_brm = {'ブロモ酢酸エチル': (1, 2), '1-ブロモヘキサン': (5, 10), 'トランス-4-ブロモ-2-ブテン酸エチル': (1, 4), 'シス-4-ブロモ-2-ブテン酸エチル': (1, 4)}
res = {}
for a, (pa, olda) in ald.items():
    for b, (pb, stab, oldb) in brm.items():
        if (pa, pb) not in valid_pairs: continue
        A, B = span_ald[a], span_brm[b]
        cover = set(range(A[0], A[1]+1)) | set(range(B[0], B[1]+1))
        if cover != set(range(1, 11)) or len(set(range(A[0], A[1]+1)) & set(range(B[0], B[1]+1))): continue
        new = valid_pairs[(pa, pb)]
        conf = {new: 'E' if stab else 'Z'}
        for o in (olda, oldb):
            if o: conf[o[0]] = o[1]
        res[(a, b)] = (conf['C2=C3'], conf['C4=C5'])
for k, v in res.items(): print(k, v)
products = set(res.values())
assert ('E', 'Z') in products and ('Z', 'Z') not in products
assert [k for k, v in res.items() if v == ('E', 'Z')] == [('シス-2-オクテナール', 'ブロモ酢酸エチル'), ('トランス-4-オキソ-2-ブテン酸エチル', '1-ブロモヘキサン')] or \
       sorted(k for k, v in res.items() if v == ('E', 'Z')) == sorted([('シス-2-オクテナール', 'ブロモ酢酸エチル'), ('トランス-4-オキソ-2-ブテン酸エチル', '1-ブロモヘキサン')])
assert [k for k, v in res.items() if v == ('Z', 'E')] == [('ヘキサナール', 'シス-4-ブロモ-2-ブテン酸エチル')]
assert res[('ヘキサナール', 'トランス-4-ブロモ-2-ブテン酸エチル')] == ('E', 'E')
assert len(res) == 5
# (5)
from sympy import Rational as R
MT = 12*12 + 20*1 + 2*16; MP = 18*12 + 15 + 16 + 31
assert (MT, MP) == (196, 278)
m = R(980, 100) / MT * MP
print('Ph3PO =', float(m), '→', sig(float(m), 2)); assert sig(float(m), 2) == 14
fr = R(MT, MT + MP) * 100
print('割合 =', float(fr)); assert round(float(fr)) == 41
print('ALL OK')
