# 検証スクリプト（グループ g10）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/POL-RUB-023/item.json'
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
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡 ハロゲン・貴ガス・水素 酸素・硫黄 窒素・リン 炭素・ケイ素 アルカリ金属・アルカリ土類金属 アルミニウム・亜鉛・スズ・鉛 遷移元素（鉄・銅・銀・クロム・マンガン） 金属イオンの沈殿・錯イオン・系統分析 気体の製法と性質 有機化合物の構造と元素分析 脂肪族炭化水素 アルコール・アルデヒド・ケトン カルボン酸・エステル・油脂 芳香族化合物 有機化合物の検出反応と分離 高分子の基礎と合成樹脂 合成繊維 ゴム・機能性高分子 糖類 アミノ酸・タンパク質・核酸'.split()
assert it['field'] in FIELDS, it['field']
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert it['solution'].rstrip().endswith('\\hfill$\\square$')

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
from collections import Counter
def ozonolysis(smi):
    """C=C をすべて切って両側を C=O にする（鎖の続き * は H とみなさず断片の外へ）"""
    m = Chem.RWMol(Chem.MolFromSmiles(smi))
    dbs = [b for b in m.GetBonds() if b.GetBondTypeAsDouble() == 2 and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C']
    pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in dbs]
    for a, b in pairs:
        m.RemoveBond(a, b)
        for x in (a, b):
            o = m.AddAtom(Chem.Atom(8)); m.AddBond(x, o, Chem.BondType.DOUBLE)
    for at in m.GetAtoms(): at.SetNoImplicit(False); at.SetNumExplicitHs(0)
    Chem.SanitizeMol(m)
    frags = Chem.GetMolFrags(m, asMols=True)
    return Counter(Chem.MolToSmiles(f) for f in frags if '*' not in Chem.MolToSmiles(f))
# 4 単位の内部断片を数える（両端の * を含む断片は除く）
a = ozonolysis('*CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CC*')
b = ozonolysis('*' + ''.join('CC(C)=CC' if i % 2 == 0 else 'CC=C(C)C' for i in range(20)) + '*')   # 20 単位，向きを交互に
c34 = ozonolysis('*CC(C(C)=C)CC(C(C)=C)*')
c12 = ozonolysis('*CC(C)(C=C)CC(C)(C=C)*')
print(a, b, c34, c12)
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
assert set(a) == {can('CC(=O)CCC=O')}
assert set(b) == {can('O=CCCC=O'), can('CC(=O)CCC(C)=O')} and abs(b[can('O=CCCC=O')] - b[can('CC(=O)CCC(C)=O')]) <= 1   # 長い鎖で 1:1（端の 1 個の差のみ）
assert set(c34) == {'C=O'} and set(c12) == {'C=O'}
for s, f in [('CC(=O)CCC=O', 'C5H8O2'), ('O=CCCC=O', 'C4H6O2'), ('CC(=O)CCC(C)=O', 'C6H10O2')]:
    assert D.CalcMolFormula(Chem.MolFromSmiles(s)) == f
# ヨードホルム（CH3-C(=O)-）・銀鏡（CHO）
io = Chem.MolFromSmarts('[CH3]C(=O)[#6,#1]'); ald = Chem.MolFromSmarts('[CX3H1](=O)')
for s, exp in [('CC(=O)CCC=O', (True, True)), ('O=CCCC=O', (False, True)), ('CC(=O)CCC(C)=O', (True, False))]:
    m = Chem.MolFromSmiles(s); r = (m.HasSubstructMatch(io), m.HasSubstructMatch(ald)); print(s, r); assert r == exp
import os
d = os.path.dirname(ITEM)
for fn in ('problem.svg', 'solution.svg'):
    t = open(os.path.join(d, fn), encoding='utf-8').read()
    assert 'viewBox' in t and '<script' not in t and 'foreignObject' not in t and 'href' not in t
print('OK')
